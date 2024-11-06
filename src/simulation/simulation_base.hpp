#pragma once

#include "simulation.hpp"
#include "spatial_grid.hpp"

void update_parameters(spatial_grid_base& grid, sph_parameters_structure const& sph_parameters);
void update_density(spatial_grid_base& grid, float h);
void update_pressure(cgp::numarray<particle_element>& particles, float rho0, float stiffness);
void update_force(spatial_grid_base& grid, sph_parameters_structure const& sph_parameters);