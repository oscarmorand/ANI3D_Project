#pragma once

#include "cgp/cgp.hpp"
#include "simulation/simulation.hpp"

// Parametric function defined as a sum of blobs-like primitives
//  f(p) = sa exp(-||p-pa||^2) + sb exp(-||p-pb||^2) + sc exp(-||p-pc||^2) + noise(p)
//   with noise: a Perlin noise
// The operator()(vec3 p) allows to query a value of the function at arbitrary point in space
struct field_function_structure {

	// Query the value of the function at any point p
	float operator()(cgp::vec3 const& p, cgp::vec3 const& particle, float radius, float radius2, float  radius22) const;


	// ***************************//
	// Parameters
	// ***************************//

	// The parameters of the Perlin noise
	float noise_magnitude   = 0.5f; // Magnitude of the noise
	float noise_offset      = 0.0f; // An offset in the parametric domain (get a different value of noise with same parameters)
	float noise_scale       = 1.0f; // Scale in the parametric domain
	int noise_octave        = 5;    // Maximum number of octave
	float noise_persistance = 0.3f; // Persistence of Perlin noise
};