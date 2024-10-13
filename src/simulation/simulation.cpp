#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Laplacian of W_viscosity
    float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
    return 45.0f/(3.14159f*std::pow(h,6)) * (h-r);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    //  Gradient of W_spiky
    vec3 const& p = p_i-p_j;
    float const r = norm(p);
    assert_cgp_no_msg(r<=h);
    return -45.0f/(3.14159f*std::pow(h,6)) * (p/r) * std::pow(h-r, 2.0f);
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}

void update_parameters(numarray<std::shared_ptr<particle_element>>& particles, grid_container& grid, float h)
{
    int const N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float m = 0.0;
        float nu = 0.0;
        vec3 color{0.0,0.0,0.0};
        float div = 0.0;

        auto fluid_i = particles[i]->fluid_type;
        auto p_i = particles[i]->p;

        auto block = particles[i]->block;
        auto neighbours = grid.get_neighbours_blocks(block);

        for(int idx=0; idx<neighbours.size(); idx++) {
            int neighbour_id = neighbours.at(idx);

            if (neighbour_id != -1) {
                auto neighbour = grid.at(neighbour_id);

                for (int j = 0; j < neighbour->size(); j++) {
                    auto const& p_j = neighbour->at(j)->p;
                    float const r = norm(p_i - p_j);
                    if (r < 2*h/3)
                    {
                        auto fluid_j = neighbour->at(j)->fluid_type;

                        if (fluid_i == fluid_j || fluid_i->soluble_classes.find(fluid_j) != fluid_i->soluble_classes.end()) {
                            m += neighbour->at(j)->m;
                            nu += neighbour->at(j)->nu;
                            color += neighbour->at(j)->color;
                            div += 1.0;
                        }
                    }
                }
            }
        }

        if (div > 0.0)
        {
            particles[i]->m = m / div;
            particles[i]->nu = nu / div;
            particles[i]->color = color / div;
        }
    }
}

void update_density(numarray<std::shared_ptr<particle_element>>& particles, grid_container& grid, float h)
{
    // Compute the density value (particles[i].rho) at each particle position

    int const N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float rho = 0;
        auto const& p_i = particles[i]->p;

        auto block = particles[i]->block;
        auto neighbours = grid.get_neighbours_blocks(block);

        for(int idx=0; idx<neighbours.size(); idx++) {
             int neighbour_id = neighbours.at(idx);

            if (neighbour_id != -1) {
                auto neighbour = grid.at(neighbour_id);

                for (int j = 0; j < neighbour->size(); j++) {
                    auto const& p_j = neighbour->at(j)->p;
                    float const m_j = neighbour->at(j)->m;

                    float const r = norm(p_i - p_j);
                    if (r<h)
                    {
                        float const w = W_density(p_i, p_j, h);
                        rho += m_j * w;
                    }
                }
            }       
        }

        particles[i]->rho = rho;
    }
}

// Convert the particle density to pressure
void update_pressure(numarray<std::shared_ptr<particle_element>>& particles, float rho0, float stiffness)
{
	const int N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
        particles[i]->pressure = density_to_pressure(particles[i]->rho, rho0, stiffness);
}

/*
// Compute the forces and update the acceleration of the particles
void update_force(numarray<std::shared_ptr<particle_element>>& particles, grid_container& grid, float h)
{
	// gravity
    const int N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        auto& particle_i = particles[i];

        float m_i = particle_i->m;
        particle_i->f = m_i * vec3{0,-9.81f,0};

        vec3 F_pressure{0,0,0};
        vec3 F_viscosity{0,0,0};

        auto p_i = particle_i->p;
        float pr_i = particle_i->pressure;
        float rho_i = particle_i->rho;
        vec3 v_i = particle_i->v;
        float nu_i = particle_i->nu;

        auto block = particle_i->block;
        auto neighbours = block->get_neighbours_blocks();

        for(int idx=0; idx<neighbours.size(); idx++) {
            int neighbour_id = neighbours.at(idx);

            if (neighbour_id != -1) {
                auto neighbour = grid.at(neighbour_id);

                for (int j = 0; j < neighbour->size(); j++) {
                    const auto& particle_j = neighbour->at(j);

                    auto p_j = particle_j->p;
                    float r = norm(p_i-p_j);
                    if (r<h && particle_i != particle_j)
                    {
                        float pr_j = particle_j->pressure;
                        float rho_j = particle_j->rho;
                        vec3 v_j = particle_j->v;
                        float m_j = particle_j->m;
                        float nu_j = particle_j->nu;

                        // Compute F_pressure
                        vec3 grad = W_gradient_pressure(p_i, p_j, h);
                        F_pressure += m_j *(pr_i+pr_j)/(2*rho_j)*grad;

                        // Compute F_viscosity
                        float lap = W_laplacian_viscosity(p_i, p_j, h);
                        F_viscosity += ((nu_j + nu_i) / 2) * m_j * (v_j - v_i)/rho_j * lap;
                    }
                }
            }
        }
        F_pressure *= -m_i / rho_i;
        F_viscosity *= m_i;

        particle_i->f += F_pressure + F_viscosity;
    }
}
*/

// Compute the forces and update the acceleration of the particles
void update_force(numarray<std::shared_ptr<particle_element>>& particles, float h)
{
	// gravity
    const int N = particles.size();

    #pragma omp parallel for
    for(int i=0; i<N; ++i)
    {
        float const m_i = particles[i]->m;
        particles[i]->f = m_i * vec3{0,-9.81f,0};

        vec3 F_pressure{0,0,0};
        vec3 F_viscosity{0,0,0};

        auto const& p_i = particles[i]->p;
        float const pr_i = particles[i]->pressure;
        float const rho_i = particles[i]->rho;
        vec3 const v_i = particles[i]->v;
        float const nu_i = particles[i]->nu;

        for(int j=0; j<N; ++j)
        {
            if (i != j) {
                auto const& p_j = particles[j]->p;
                float const pr_j = particles[j]->pressure;
                float const rho_j = particles[j]->rho;
                vec3 const v_j = particles[j]->v;
                float const m_j = particles[j]->m;
                float const nu_j = particles[j]->nu;

                float const r = norm(p_i-p_j);
                if (r<h)
                {
                    // Compute F_pressure
                    vec3 const grad = W_gradient_pressure(p_i, p_j, h);
                    F_pressure += m_j *(pr_i+pr_j)/(2*rho_j)*grad;

                    // Compute F_viscosity
                    float const lap = W_laplacian_viscosity(p_i, p_j, h);
                    F_viscosity += ((nu_j + nu_i) / 2) * m_j * (v_j - v_i)/rho_j * lap;
                }
            }
        }
        F_pressure *= -m_i / rho_i;
        F_viscosity *= m_i;

        particles[i]->f += F_pressure + F_viscosity;
    }
}

void construct_grid(numarray<std::shared_ptr<particle_element>>& particles, grid_container& grid)
{
    grid.clear();

    int const N = particles.size();

    //#pragma omp parallel for FIXME
    for(int i=0; i<N; ++i)
    {
		grid.add_particle(particles[i]);
    }
}

void simulate(float dt, numarray<std::shared_ptr<particle_element>>& particles, grid_container& grid, sph_parameters_structure const& sph_parameters)
{
    construct_grid(particles, grid);

    update_parameters(particles, grid, sph_parameters.h);
    update_density(particles, grid, sph_parameters.h);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    //update_force(particles, grid, sph_parameters.h);  // Update forces
    update_force(particles, sph_parameters.h);

	// Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();

    // Collision
    float const epsilon = 1e-3f;

    #pragma omp parallel for
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k]->p;
		vec3& v = particles[k]->v;
		vec3& f = particles[k]->f;
        float const m = particles[k]->m;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.y>1 )  {p.y =  1-epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.z<-0.01 ) {p.z = -0.01+epsilon*rand_uniform();  v.z *= -0.5f;}
        if( p.z>0.01 )  {p.z =  0.01-epsilon*rand_uniform();  v.z *= -0.5f;}
	}
}