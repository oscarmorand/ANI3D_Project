#pragma once

#include "cgp/cgp.hpp"

#include <unordered_set>
#include <memory>

struct fluid_class
{
    std::unordered_set<std::shared_ptr<fluid_class>> soluble_classes;

    std::string fluid_name;

    float base_m;
    float base_nu;
    cgp::vec3 base_color;

    fluid_class() : fluid_name(""), base_m(0), base_nu(0), base_color{0,0,0} {};

    fluid_class(std::string name, float m, float nu, cgp::vec3 color) : fluid_name(name), base_m(m), base_nu(nu), base_color(color) {};
};

// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    cgp::vec3 external_forces; 

    float m;        // mass of the particle
    float rho;      // density at this particle position
    float pressure; // pressure at this particle position
    float nu;       // viscosity

    std::shared_ptr<fluid_class> fluid_type;

    cgp::vec3 color; // color of the particle

    int cell_index; // Index of the cell in the spatial grid

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
    
    float gravity_strength = 9.81f;

    float fluid_mixing_rate = 1.0f;
};