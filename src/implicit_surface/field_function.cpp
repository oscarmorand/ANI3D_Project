#include "implicit_surface/field_function.hpp"

using namespace cgp;

// Parameterization Gaussian centered at point p0
static float gaussian(vec3 const& p, vec3 const& p0, float sigma)
{
	float const d = norm(p - p0);
	return d < sigma ? 1.0f / (d + 1e-6f) : 0.0f;
}

static float quartic_falloff(cgp::vec3 const& p, cgp::vec3 const& p0, float radius)
{
    float const d = norm(p - p0);
    if (d >= radius) return 0.0f;
    float t = d / radius;
    return (1.0f - t) * (1.0f - t) * (1.0f - t) * (1.0f - t);  // Quartic falloff
}

static float laplacian_of_gaussian(cgp::vec3 const& p, cgp::vec3 const& p0, float sigma2, float sigma22)
{
    float const d = norm(p - p0);
    return (1 - (d * d) / (sigma2)) * std::exp(-(d * d) / (2 * sigma2));
}

static float exponential_decay(cgp::vec3 const& p, cgp::vec3 const& p0, float sigma)
{
    float const d = norm(p - p0);
    return std::exp(-d / sigma);
}

static float cubic_falloff(cgp::vec3 const& p, cgp::vec3 const& p0, float radius)
{
    float const d = norm(p - p0);
    if (d >= radius) return 0.0f;
    float t = d / radius;
    return (1.0f - t) * (1.0f - t) * (1.0f - t);  // Cubic falloff
}

float field_function_structure::operator()(cgp::vec3 const& cell, cgp::vec3 const& particle, float radius2, float radius22) const
{
	float value = 0.0f;
	value += laplacian_of_gaussian(particle, cell, radius2, radius22);
	return value;
}