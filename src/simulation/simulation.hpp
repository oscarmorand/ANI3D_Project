#pragma once

#include "cgp/cgp.hpp"

#include "particle_element.hpp"
#include "fluid_class.hpp"
#include "grid_container.hpp"

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;  
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};

void simulate(float dt, cgp::numarray<std::shared_ptr<particle_element>>& particles, grid_container& particles_grid, sph_parameters_structure const& sph_parameters);