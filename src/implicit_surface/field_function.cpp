#include "implicit_surface/field_function.hpp"

using namespace cgp;

// Parameterization Gaussian centered at point p0
static float gaussian(vec3 const& p, vec3 const& p0, float sigma)
{
	float const d = norm(p - p0);
	float const value = std::exp(-(d * d) / (sigma * sigma));
	return value;
}

float field_function_structure::operator()(cgp::vec3 const& p, cgp::numarray<particle_element> particules) const
{
	float value = 0.0f;
	#pragma omp parallel for
	for (int i = 0; i < particules.size(); i++) {
		vec3 const& p0 = particules[i].p;
		float const d = norm(p - p0);
		float const v = std::exp(-(d * d) / (0.01));
		value += v;
	}

	return value;
}