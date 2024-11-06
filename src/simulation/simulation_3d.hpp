#pragma once

#include "simulation.hpp"
#include "spatial_grid.hpp"

void simulate_3d(float dt, cgp::numarray<particle_element>& particles, spatial_grid_base& grid, sph_parameters_structure const& sph_parameters);