#pragma once

#include "simulation.hpp"
#include "spatial_grid.hpp"

void simulate_2d(float dt, cgp::numarray<particle_element>& particles, spatial_grid& grid, sph_parameters_structure const& sph_parameters);