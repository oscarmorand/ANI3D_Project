#include "scene.hpp"

void scene_structure::new_fluid_gui()
{
    if (!gui.new_fluid_creation) {
        ImGui::SameLine();
        bool new_fluid_pressed = ImGui::Button("New fluid");
        if (new_fluid_pressed) {
            gui.new_fluid_creation = true;

			gui.new_fluid_name[0] = fluid_classes.size() + '0';
			gui.new_fluid_name[1] = '\0';
			gui.new_fluid_color = { 1.0, 1.0, 1.0 };
			gui.new_fluid_mass = 1.0f;
			gui.new_fluid_viscosity = 0.01f;

			/*
			gui.new_fluid_soluble.clear();
            gui.new_fluid_soluble.resize(fluid_classes.size());

            for (int i = 0; i < fluid_classes.size(); i++) {
                gui.new_fluid_soluble[i] = std::make_shared<bool>(false);
            }
			*/
        }
    }
    if (gui.new_fluid_creation)
    {
        ImGui::InputText("Fluid name", gui.new_fluid_name, 64);
        ImGui::SliderFloat("Mass", &gui.new_fluid_mass, 0.01f, 3.0f, "%0.2f");
        ImGui::SliderFloat("Viscosity", &gui.new_fluid_viscosity, 0.0f, 0.1f, "%0.3f");
        ImGui::ColorPicker3("Color", &gui.new_fluid_color[0]);

		/*
		for (int i = 0; i < fluid_classes.size(); i++) {
			ImGui::Checkbox(fluid_classes[i]->fluid_name.c_str(), gui.new_fluid_soluble[i].get());
		}
		*/
		
        bool create_pressed = ImGui::Button("Create fluid"); ImGui::SameLine();
        bool cancel_pressed = ImGui::Button("Cancel");

        if (create_pressed) {
            gui.new_fluid_creation = false;

            auto new_fluid = std::make_shared<fluid_class>(gui.new_fluid_name, gui.new_fluid_mass, gui.new_fluid_viscosity, gui.new_fluid_color);

			/*
			for (int i = 0; i < fluid_classes.size(); i++) {
				if (*(gui.new_fluid_soluble[i])) {
                    new_fluid->soluble_classes.insert(fluid_classes[i]);
                    fluid_classes[i]->soluble_classes.insert(new_fluid);
                }
			}
			*/

			fluid_classes.push_back(new_fluid);
        }

        if (cancel_pressed) {
            gui.new_fluid_creation = false;
        }

        ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
    }
}

void scene_structure::display_gui()
{
	bool change_dimension = ImGui::RadioButton("2D", &dimension, DIM_2D); ImGui::SameLine();
	change_dimension |= ImGui::RadioButton("3D", &dimension, DIM_3D);
	if (change_dimension)
	{
		initialize_sph();
		camera_control = camera_controller();
		camera_control.initialize(inputs, window);
		if (dimension == DIM_2D)
		{
			gui.display_color = true;
			cam_projection = cgp::camera_projection::build_orthographic(-1.1f, 1.1f, -1.1f, 1.1f, -10, 10);
		}
		else
		{
			gui.display_color = false;
			cam_projection = cgp::camera_projection::build_perspective(50.0f * 3.14f / 180, 1.0f, 0.1f, 100.0f);
		}
	}

	bool const play_pause = ImGui::Button(timer.is_running() ? "Pause" : "Play");
	if (play_pause)
	{
		if (timer.is_running())
			timer.stop();
		else
			timer.start();
	}

	bool const restart = ImGui::Button("Restart");
	if (restart)
	{
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
	if (gravity_preset)
	{
		sph_parameters.gravity_strength = 9.81f;
	}
	gravity_preset = ImGui::Button("Moon"); ImGui::SameLine();
	if (gravity_preset)
	{
		sph_parameters.gravity_strength = 1.62f;
	}
	gravity_preset = ImGui::Button("Jupiter");
	if (gravity_preset)
	{
		sph_parameters.gravity_strength = 24.79f;
	}
	ImGui::SliderFloat("Neighbour radius", &sph_parameters.h, 0.1f, 0.3f, "%0.3f");
	ImGui::SliderFloat("Fluid mixing rate", &sph_parameters.fluid_mixing_rate, 0.0f, 1.0f, "%0.2f");

	ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();
	ImGui::Text("Right-click action");
	ImGui::RadioButton("Spawn particles", &gui.right_click_action, SPAWN_PARTICLES); ImGui::SameLine();
	ImGui::RadioButton("Remove particles", &gui.right_click_action, REMOVE_PARTICLES); ImGui::SameLine();
	ImGui::RadioButton("Add repulsive force", &gui.right_click_action, ADD_RADIAL_FORCE);
	ImGui::RadioButton("Add vortex", &gui.right_click_action, ADD_VORTEX_FORCE); ImGui::SameLine();
	ImGui::RadioButton("Add attractive force", &gui.right_click_action, ADD_GRAVITY_FORCE);

	if (gui.right_click_action == SPAWN_PARTICLES)
	{
		int i = 0;
		for (auto const& fluid : fluid_classes)
		{
			ImGui::RadioButton(fluid->fluid_name.c_str(), &gui.spawn_particle_type, i);
			if (i < fluid_classes.size() - 1)
				ImGui::SameLine();
			i += 1;
		}
        new_fluid_gui();

		ImGui::SliderInt("Number of particles", &gui.spawn_particle_number, 1, 1000);
		if (dimension == DIM_2D)
			ImGui::SliderFloat("Radius of spawn disk", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
		else
			ImGui::SliderFloat("Radius of spawn sphere", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
	}
	else if (gui.right_click_action == REMOVE_PARTICLES)
	{
		if (dimension == DIM_2D)
			ImGui::SliderFloat("Radius of deletion disk", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
		else
			ImGui::SliderFloat("Radius of deletion sphere", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
	}
	else if (gui.right_click_action == ADD_RADIAL_FORCE || gui.right_click_action == ADD_VORTEX_FORCE || gui.right_click_action == ADD_GRAVITY_FORCE)
	{
		if (dimension == DIM_2D)
			ImGui::SliderFloat("Radius of force disk", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
		else
			ImGui::SliderFloat("Radius of force sphere", &gui.spawn_particle_radius, 0.01f, 1.0f, "%0.2f");
		ImGui::SliderFloat("Force strength", &gui.force_strength, 0.01f, 3.0f, "%0.2f");
		if (gui.right_click_action == ADD_VORTEX_FORCE)
		{
			ImGui::RadioButton("Clockwise", &gui.vortex_direction, CLOCKWISE); ImGui::SameLine();
			ImGui::RadioButton("Counter-clockwise", &gui.vortex_direction, COUNTER_CLOCKWISE);
		}
	}

	ImGui::Spacing();ImGui::Spacing();ImGui::Spacing();
	ImGui::Text("Display options");
	if (dimension == DIM_2D) {
		ImGui::Checkbox("Color", &gui.display_color);

		if (gui.display_color)
		{
			bool grid_size_changed = ImGui::SliderInt("Grid size", &gui.grid_size, 1, 100);
			if (grid_size_changed)
			{
				field.resize(gui.grid_size, gui.grid_size);
				field_quad.texture.initialize_texture_2d_on_gpu(field);
			}
			ImGui::RadioButton("Fluid color", &gui.color_type, FLUID_COLOR);ImGui::SameLine();
			ImGui::RadioButton("Velocity", &gui.color_type, VELOCITY);ImGui::SameLine();
			ImGui::RadioButton("Density", &gui.color_type, DENSITY);

			if (gui.color_type == VELOCITY || gui.color_type == DENSITY)
			{
				ImGui::SliderFloat("Smoothing radius", &gui.color_smoothing_radius, 0.0f, 0.3f);

				bool min_changed = false;
				bool max_changed = false;
				if (gui.color_type == VELOCITY)
				{
					min_changed = ImGui::SliderFloat("Velocity min", &gui.threshold_min, 0.0f, 1.0f);
					max_changed = ImGui::SliderFloat("Velocity max", &gui.threshold_max, 0.0f, 1.0f);
				}
				else if (gui.color_type == DENSITY)
				{
					min_changed = ImGui::SliderFloat("Density min", &gui.threshold_min, 0.0f, 1.0f);
					max_changed = ImGui::SliderFloat("Density max", &gui.threshold_max, 0.0f, 1.0f);
				}

				if (gui.threshold_min > gui.threshold_max) {
					if (min_changed)
						gui.threshold_max = gui.threshold_min + 0.01f;
					if (max_changed)
						gui.threshold_min = gui.threshold_max - 0.01f;
				}

				ImGui::ColorPicker3("Color min", &gui.color_min[0]);ImGui::SameLine();
				ImGui::ColorPicker3("Color max", &gui.color_max[0]);
			}
		}
	}
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Radius", &gui.display_radius);
}