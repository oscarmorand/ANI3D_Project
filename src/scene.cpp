#include "scene.hpp"

using namespace cgp;


void scene_structure::initialize()
{
	if (dimension == DIM_2D)
		cam_projection = cgp::camera_projection::build_orthographic(-1.1f, 1.1f, -1.1f, 1.1f, -10, 10);
	else
		cam_projection = cgp::camera_projection::build_perspective(50.0f * 3.14f / 180, 1.0f, 0.1f, 100.0f);
	cam_projection.update_aspect_ratio(window.aspect_ratio());
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.look_at({ 0.0f, 0.0f, 2.0f }, {0,0,0}, {0,1,0});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	field.resize(50, 50);
	field_quad.initialize_data_on_gpu(mesh_primitive_quadrangle({ -1,-1,0 }, { 1,-1,0 }, { 1,1,0 }, { -1,1,0 }) );
	field_quad.material.phong = { 1,0,0 };
	field_quad.texture.initialize_texture_2d_on_gpu(field);

	initialize_sph();
	sphere_particle.initialize_data_on_gpu(mesh_primitive_sphere(1.0,{0,0,0},10,10));
	sphere_particle.model.scaling = 0.01f;
	curve_visual.color = { 1,0,0 };
	curve_visual.initialize_data_on_gpu(curve_primitive_circle());
}

void scene_structure::spawn_particle(vec3 const& pos, fluid_type_enum fluid_type)
{
	particle_element particle;
	particle.p = pos;
	particle.v = { 0,0,0 };
	particle.f = { 0,0,0 };
	
	switch (fluid_type) {
		case WATER:
			particle.m = 1.0f;
			particle.nu = 0.0005f;
			particle.color = { 0.1, 0.3, 1.0 };
			break;
		case MILK:
			particle.m = 1.0f;
			particle.nu = 0.0005f;
			particle.color = { 1.0, 1.0, 1.0 };
			break;
		case OIL:
			particle.m = 0.1f;
			particle.nu = 0.05f;
			particle.color = { 1.0, 1.0, 0.3 };
			break;
	}
	particle.fluid_type = fluid_classes[fluid_type];

	particles.push_back(particle);
}

void scene_structure::spawn_particles_in_disk(vec3 const& center, float radius, int N, fluid_type_enum fluid_type)
{
    for (int k = 0; k < N; ++k) {
        vec2 const u = { rand_uniform(), rand_uniform() };
        vec3 const p = { center.x + radius * (2 * u.x - 1), center.y + radius * (2 * u.y - 1), 0};
        spawn_particle(p, fluid_type);
    }
}

void scene_structure::spawn_particles_in_sphere(vec3 const& center, float radius, int N, fluid_type_enum fluid_type)
{
    for (int k = 0; k < N; ++k) {
        vec3 const u = { rand_uniform(), rand_uniform(), rand_uniform() };
        vec3 const p = { center.x + radius * (2 * u.x - 1), center.y + radius * (2 * u.y - 1), center.z + radius * (2 * u.z - 1) };
        spawn_particle(p, fluid_type);
    }
}

void scene_structure::initialize_sph()
{
	// Initial particle spacing (relative to h)
	float const c = 1.0f;
	float const h = sph_parameters.h;

	auto water = std::make_shared<fluid_class>();
	auto milk = std::make_shared<fluid_class>();
	auto oil = std::make_shared<fluid_class>();

	water->soluble_classes.insert(milk);
	milk->soluble_classes.insert(water);

	fluid_classes[WATER] = water;
	fluid_classes[MILK] = milk;
	fluid_classes[OIL] = oil;

	// Fill a square with particles
	particles.clear();
	for (float x = h; x < 1.0f - h; x = x + c * h)
	{
		for (float y = -1.0f + h; y < 1.0f - h; y = y + c * h)
		{
			for (float z = -1.0f + h; z < 1.0f - h; z = z + c * h)
			{
				vec3 pos = {0, 0, 0};
				if (dimension == DIM_2D)
					pos = { x + h / 8.0 * rand_uniform(), y + h / 8.0 * rand_uniform(), 0 };
				else
					pos = { x + h / 8.0 * rand_uniform(), y + h / 8.0 * rand_uniform(), z + h / 8.0 * rand_uniform() };

				fluid_type_enum fluid_type = WATER;
				float m = rand_uniform();
				if (m < 0.3f) {
					fluid_type = WATER;
				}
				else if (m < 0.6f) {
					fluid_type = OIL;
				}
				else {
					fluid_type = MILK;
				}

				spawn_particle(pos, fluid_type);
			}
		}
	}
}

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	if (timer.is_running()) {
		timer.update(); // update the timer to the current elapsed time
		float const dt = 0.005f * timer.scale;

		if (dimension == DIM_2D)
			simulate_2d(dt, particles, sph_parameters);
		else
			simulate_3d(dt, particles, sph_parameters);
	}

	if (gui.display_particles) {
		for (int k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.model.translation = p;
			sphere_particle.material.color = particles[k].color;
			draw(sphere_particle, environment);
		}
	}

	if (gui.display_radius) {
		curve_visual.model.scaling = sph_parameters.h;
		for (int k = 0; k < particles.size(); k += 10) {
			curve_visual.model.translation = particles[k].p;
			draw(curve_visual, environment);
		}
	}

	if (gui.display_color) {
		update_field_color();
		field_quad.texture.update(field);
		draw(field_quad, environment);
	}
}

void scene_structure::display_gui()
{
	bool change_dimension = ImGui::RadioButton("2D", &dimension, DIM_2D); ImGui::SameLine();
	change_dimension |= ImGui::RadioButton("3D", &dimension, DIM_3D);
	if (change_dimension) {
		initialize_sph();
		if (dimension == DIM_2D) {
			gui.display_color = true;
			cam_projection = cgp::camera_projection::build_orthographic(-1.1f, 1.1f, -1.1f, 1.1f, -10, 10);
		}
		else {
			gui.display_color = false;
			cam_projection = cgp::camera_projection::build_perspective(50.0f * 3.14f / 180, 1.0f, 0.1f, 100.0f);
		}
	}

	bool const play_pause = ImGui::Button(timer.is_running() ? "Pause" : "Play");
	if (play_pause) {
		if (timer.is_running())
			timer.stop();
		else
			timer.start();
	}

	bool const restart = ImGui::Button("Restart");
	if (restart) {
		initialize_sph();
	}
	
	bool const clear_particles = ImGui::Button("Clear particles");
	if (clear_particles)
		particles.clear();

	ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
	ImGui::Text("Simulation parameters");
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");
	ImGui::SliderFloat("Gravity", &sph_parameters.gravity_strength, 0.0f, 30.0f, "%0.1f");
	ImGui::Text("Gravity Preset"); ImGui::SameLine();
	bool gravity_preset = ImGui::Button("Earth"); ImGui::SameLine();
	if (gravity_preset) {
		sph_parameters.gravity_strength = 9.81f;
	}
	gravity_preset = ImGui::Button("Moon"); ImGui::SameLine();
	if (gravity_preset) {
		sph_parameters.gravity_strength = 1.62f;
	}
	gravity_preset = ImGui::Button("Jupiter");
	if (gravity_preset) {
		sph_parameters.gravity_strength = 24.79f;
	}
	ImGui::SliderFloat("Neighbour radius", &sph_parameters.h, 0.1f, 0.3f, "%0.3f");
	ImGui::SliderFloat("Fluid mixing rate", &sph_parameters.fluid_mixing_rate, 0.0f, 1.0f, "%0.2f");

	ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
	ImGui::Text("Right-click action");
	ImGui::RadioButton("Spawn particles", &gui.right_click_action, SPAWN_PARTICLES); ImGui::SameLine();
	ImGui::RadioButton("Add force", &gui.right_click_action, ADD_FORCE);

	if (gui.right_click_action == SPAWN_PARTICLES) {
		ImGui::RadioButton("Water", &gui.spawn_particle_type, WATER); ImGui::SameLine();
		ImGui::RadioButton("Milk", &gui.spawn_particle_type, MILK); ImGui::SameLine();
		ImGui::RadioButton("Oil", &gui.spawn_particle_type, OIL);

		ImGui::SliderInt("Number of particles", &gui.spawn_particle_number, 1, 1000);
		if (dimension == DIM_2D)
			ImGui::SliderFloat("Radius of spawn disk", &gui.spawn_particle_radius, 0.01f, 0.5f, "%0.2f");
		else
			ImGui::SliderFloat("Radius of spawn sphere", &gui.spawn_particle_radius, 0.01f, 0.5f, "%0.2f");
	}
	else {
		// TODO
	}

	ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
	ImGui::Text("Display options");
	if (dimension == DIM_2D)
		ImGui::Checkbox("Color", &gui.display_color);
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Radius", &gui.display_radius);

	if (gui.display_color) {
		ImGui::Text("2D field color");
		ImGui::RadioButton("Fluid color", &gui.color_type, FLUID_COLOR); ImGui::SameLine();
		ImGui::RadioButton("Velocity", &gui.color_type, VELOCITY);

		if (gui.color_type == VELOCITY) {
		ImGui::SliderFloat("Velocity min", &gui.threshold_min, 0.0f, 10.0f);
		ImGui::SliderFloat("Velocity max", &gui.threshold_max, 0.0f, 10.0f);

		ImGui::ColorPicker3("Color min", &gui.color_min[0]); ImGui::SameLine();
		ImGui::ColorPicker3("Color max", &gui.color_max[0]);
		}
	}
}

vec3 color_interpolation(vec3 const& c0, vec3 const& c1, float x0, float x1, float x)
{
	float const alpha = (x - x0) / (x1 - x0);
	return (1 - alpha) * c0 + alpha * c1;
}

vec3 color_clamp(vec3 const& v, float vmin, float vmax)
{
	return { clamp(v.x, vmin, vmax), clamp(v.y, vmin, vmax), clamp(v.z, vmin, vmax) };
}

vec3 scene_structure::get_particle_color(particle_element const& particle)
{
	vec3 res = { 1,1,1 };

	switch (gui.color_type) {
		case FLUID_COLOR: {
			res = particle.color;
			break;
		}
		case VELOCITY: {
			vec3 const v = particle.v;
			res = color_clamp(
				color_interpolation(
					gui.color_min, gui.color_max, gui.threshold_min, gui.threshold_max, norm(v)
				), 0.0, 1.0
			);
			break;
		}
	}
	return res;
}

void scene_structure::update_field_color()
{
	field.fill({ 1,1,1 });
	float const d = 0.1f;
	int const Nf = int(field.dimension.x);

	#pragma omp parallel for
	for (int kx = 0; kx < Nf; ++kx) {
		for (int ky = 0; ky < Nf; ++ky) {

			vec3 full_color = { 0, 0, 0 };
			vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };

			float min_dist = d;

			for (size_t k = 0; k < particles.size(); ++k) {
				vec3 const& color = get_particle_color(particles[k]);
				vec3 const& pi = particles[k].p;
				float const r = norm(pi - p0);

				if (r < min_dist) {
					min_dist = r;
					full_color = color;
				}
			}

			field(kx, Nf - 1 - ky) = full_color;
		}
	}
}

void scene_structure::mouse_move_event()
{
	// Do not mode the camera
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
}

void scene_structure::mouse_click_event()
{
	if (inputs.mouse.click.left) { // Movements
		camera_control.action_mouse_click(environment.camera_view);
	}
	if (inputs.mouse.click.right) { // Special action
		if (gui.right_click_action == SPAWN_PARTICLES) {
			if (dimension == DIM_2D) {
				vec2 const cursor = inputs.mouse.position.current;
				vec3 const p = { cursor.x, cursor.y, 0 };
				spawn_particles_in_disk(p, gui.spawn_particle_radius, gui.spawn_particle_number, fluid_type_enum(gui.spawn_particle_type));
			}
			else {
				// TODO
			}
		}
	}
}

void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

