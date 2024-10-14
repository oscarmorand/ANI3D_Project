#include "simulation_2d.hpp"

using namespace cgp;

// Convert a density value to a pressure
static float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

static float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Laplacian of W_viscosity
    float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
    return 45.0f/(3.14159f*std::pow(h,6)) * (h-r);
}

static vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Gradient of W_spiky
    vec3 const& p = p_i-p_j;
    float const r = norm(p);
    assert_cgp_no_msg(r<=h);
    return -45.0f/(3.14159f*std::pow(h,6)) * (p/r) * std::pow(h-r, 2.0f);
}

static float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}

static void update_parameters(numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{
    int const N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float m = 0.0;
        float nu = 0.0;
        vec3 color{0.0,0.0,0.0};
        float div = 0.0;

        auto fluid_i = particles[i].fluid_type;
        const auto& p_i = particles[i].p;
        auto v_i = norm(particles[i].v);

        #pragma omp parallel for
        for(int j=0; j<N; ++j)
        {
            if (i == j) continue;

            auto const& p_j = particles[j].p;
            float const r = norm(p_i - p_j);
            if (r < sph_parameters.h)
            {
                auto fluid_j = particles[j].fluid_type;

                if (fluid_i == fluid_j || fluid_i->soluble_classes.find(fluid_j) != fluid_i->soluble_classes.end()) {
                    m += particles[j].m;
                    nu += particles[j].nu;
                    color += particles[j].color;
                    div += 1.0;
                }
            }
        }

        if (div > 0.0)
        {
            float fixed_rate = 0.1f;

            float diff_m = (m / div) - particles[i].m;
            particles[i].m += sph_parameters.fluid_mixing_rate * diff_m * fixed_rate * v_i;

            float diff_nu = (nu / div) - particles[i].nu;
            particles[i].nu += sph_parameters.fluid_mixing_rate * diff_nu * fixed_rate * v_i;

            vec3 diff_color = (color / div) - particles[i].color;
            particles[i].color += sph_parameters.fluid_mixing_rate * diff_color * fixed_rate * v_i;
        }
    }
}

static void update_density(numarray<particle_element>& particles, float h)
{
    // Compute the density value (particles[i].rho) at each particle position

    int const N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float rho = 0;
        auto const& p_i = particles[i].p;

        #pragma omp parallel for
        for(int j=0; j<N; ++j)
        {
            auto const& p_j = particles[j].p;
            float const m_j = particles[j].m;

            float const r = norm(p_i - p_j);
            if (r<h)
            {
                float const w = W_density(p_i, p_j, h);
                rho += m_j * w;
            }
        }

        particles[i].rho = rho;
    }
}

// Convert the particle density to pressure
static void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
static void update_force(numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{
    int N = particles.size();
    float h = sph_parameters.h;

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float m_i = particles[i].m;
        particles[i].f = m_i * vec3{0,-1.0f,0} * sph_parameters.gravity_strength;
    
        vec3 F_pressure{0,0,0};
        vec3 F_viscosity{0,0,0};

        vec3 const& p_i = particles[i].p;
        float pr_i = particles[i].pressure;
        float rho_i = particles[i].rho;
        vec3 const& v_i = particles[i].v;
        float nu_i = particles[i].nu;

        #pragma omp parallel for
        for(int j=0; j<N; ++j)
        {
            if (i != j) {
                vec3 const& p_j = particles[j].p;
                float pr_j = particles[j].pressure;
                float rho_j = particles[j].rho;
                vec3 const& v_j = particles[j].v;
                float m_j = particles[j].m;
                float nu_j = particles[j].nu;

                float const r = norm(p_i-p_j);
                if (r<h)
                {
                    // Compute F_pressure
                    vec3 grad = W_gradient_pressure(p_i, p_j, h);
                    F_pressure += m_j *(pr_i+pr_j)/(2*rho_j)*grad;

                    // Compute F_viscosity
                    float lap = W_laplacian_viscosity(p_i, p_j, h);
                    F_viscosity += ((nu_j + nu_i) / 2) * m_j * (v_j - v_i)/rho_j * lap;
                }
            }
        }
        F_pressure *= -m_i / rho_i;
        F_viscosity *= m_i;

        particles[i].f += F_pressure + F_viscosity;
    }
}

void simulate_2d(float dt, numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{
	// Update values
    if (sph_parameters.fluid_mixing_rate > 0.0f)
        update_parameters(particles, sph_parameters);
    update_density(particles, sph_parameters.h);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();

    // Collision
    float const epsilon = 1e-3f;

    #pragma omp parallel for
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;
        float const m = particles[k].m;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.z != 0) {p.z = 0;  v.z = 0;}
	}
}