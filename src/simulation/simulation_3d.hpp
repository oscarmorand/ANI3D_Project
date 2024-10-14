#pragma once

#include "simulation.hpp"

void simulate_3d(float dt, cgp::numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters);