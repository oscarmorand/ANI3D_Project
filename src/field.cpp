#include "scene.hpp"

vec3 color_interpolation(vec3 const &c0, vec3 const &c1, float x0, float x1, float x)
{
	float const alpha = (x - x0) / (x1 - x0);
	return (1 - alpha) * c0 + alpha * c1;
}

vec3 color_clamp(vec3 const &v, float vmin, float vmax)
{
	return {clamp(v.x, vmin, vmax), clamp(v.y, vmin, vmax), clamp(v.z, vmin, vmax)};
}

vec3 scene_structure::get_particle_color(particle_element const &particle)
{
	vec3 res = {1, 1, 1};

	switch (gui.color_type)
	{
	case FLUID_COLOR:
	{
		res = particle.color;
		break;
	}
	case VELOCITY:
	{
		vec3 const v = particle.v;
		res = color_clamp(
			color_interpolation(
				gui.color_min, gui.color_max, gui.threshold_min, gui.threshold_max, norm(v)),
			0.0, 1.0);
		break;
	}
	}
	return res;
}

void scene_structure::update_field_closest(int Nf)
{
	// Iterate over each grid cell in the field
	#pragma omp parallel for
	for (int kx = 0; kx < Nf; ++kx)
	{
		for (int ky = 0; ky < Nf; ++ky)
		{
			vec3 full_color = {0, 0, 0};  // Store the resulting color
			vec3 const p0 = {2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f};
			int cell_index = grid.compute_cell_index(p0);
			float min_dist = sph_parameters.h;

			std::vector<particle_element*> neighbors = grid.get_neighbors(cell_index);

			for (particle_element* particle : neighbors)
			{
				vec3 const& color = get_particle_color(*particle);
				vec3 const& pi = particle->p;
				float const r = norm(pi - p0);

				if (r < min_dist)
				{
					min_dist = r;
					full_color = color;
				}
			}

			field(kx, Nf - 1 - ky) = full_color;
		}
	}
}

void scene_structure::update_field_mean(int Nf)
{
	#pragma omp parallel for
	for (int kx = 0; kx < Nf; ++kx)
	{
		for (int ky = 0; ky < Nf; ++ky)
		{
			vec3 full_color = {0, 0, 0};
			int nb = 0;
			vec3 const p0 = {2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f};

			for (size_t k = 0; k < particles.size(); ++k)
			{
				vec3 const &color = get_particle_color(particles[k]);
				vec3 const &pi = particles[k].p;
				float const r = norm(pi - p0);

				if (r < gui.color_smoothing_radius) {
					full_color += color;
					nb += 1;
				}
			}

			if (nb > 0)
				full_color /= nb;
			field(kx, Nf - 1 - ky) = full_color;
		}
	}
}

void scene_structure::update_field_color()
{
	field.fill({1, 1, 1});
	int const Nf = int(field.dimension.x);

	if (gui.color_type == FLUID_COLOR)
	{
		update_field_closest(Nf);
	}
	else if (gui.color_type == VELOCITY)
	{
		update_field_mean(Nf);
	}
}