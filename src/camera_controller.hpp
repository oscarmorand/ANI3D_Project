#pragma once

#include "camera.hpp"
#include "utils.hpp"

// Specialized camera controller representing an "orbit" camera motion (= camera that rotates around a central point)
// 	- camera_controller_orbit relies on an internal camera using a quaternion to describe the camera orientation.
//  This camera controller allows to rotate the camera in arbitrary direction in 3D and models by default an ArcBall system

using namespace cgp;

struct camera_controller
{
    camera camera_model;

    // Allow to activate/deactivate the camera
	bool is_active = true;

    // Pointers to the global state of the inputs (keyboard, mouse, etc)
	input_devices* inputs = nullptr;
	// Pointer to the global state of the window
	window_structure* window = nullptr;

    void initialize(input_devices& inputs, window_structure& window);

    void action_mouse_move(mat4 &camera_matrix_view, int dim, camera_projection& cam_projection);
    void action_keyboard(mat4&) {};
    void idle_frame(mat4 &camera_matrix_view);

    void look_at(vec3 const &eye, vec3 const &center, vec3 const &up);
    void update(mat4 &camera_matrix_view);

    std::string doc_usage() const;
};
