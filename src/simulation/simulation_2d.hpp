#pragma once

#include "simulation.hpp"

void simulate_2d(float dt, cgp::numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters);