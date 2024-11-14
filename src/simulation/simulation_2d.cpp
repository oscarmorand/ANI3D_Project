#include "simulation_2d.hpp"

#include "simulation_base.hpp"

using namespace cgp;

void simulate_2d(float dt, numarray<particle_element>& particles, spatial_grid_base& grid, sph_parameters_structure const& sph_parameters)
{
	// Update values
    if (sph_parameters.fluid_mixing_rate > 0.0f)
        update_parameters(grid, sph_parameters);
    update_density(grid, sph_parameters.h);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(grid, sph_parameters);  // Update forces

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
        if( p.y>1 )  {p.y =  1-epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.z != 0) {p.z = 0;  v.z = 0;}
	}
}