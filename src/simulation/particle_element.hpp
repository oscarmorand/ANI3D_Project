#pragma once

#include "cgp/cgp.hpp"

#include "fluid_class.hpp"

#include <memory>

class block_container;

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

    std::shared_ptr<block_container> block;

    cgp::vec3 color; // color of the particle

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},m(0),rho(0),pressure(0) {}
};