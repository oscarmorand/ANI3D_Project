#pragma once

#include "cgp/cgp.hpp"

#include <unordered_set>
#include <memory>

struct fluid_class
{
    std::unordered_set<std::shared_ptr<fluid_class>> soluble_classes;
};

// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float m;        // mass of the particle
    float rho;      // density at this particle position
    float pressure; // pressure at this particle position
    float nu;       // viscosity

    std::shared_ptr<fluid_class> fluid_type;

    cgp::vec3 color; // color of the particle

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},m(0),rho(0),pressure(0) {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;  
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
    float gravity_strength = 1.0f;

    float fluid_mixing_rate = 1.0f;
};