#include "cgp/01_base/base.hpp"
#include "camera/camera_controller.hpp"

void camera_controller::initialize(input_devices& inputs_param, window_structure& window_param)
{
	inputs = &inputs_param;
	window = &window_param;
}

void camera_controller::update(mat4 &camera_matrix_view)
{
    camera_matrix_view = camera_model.matrix_view();
}

void camera_controller::action_mouse_move(mat4 &camera_matrix_view, int dim, camera_projection& cam_projection)
{
    // Preconditions
    assert_cgp_no_msg(inputs != nullptr);
    assert_cgp_no_msg(window != nullptr);
    if (!is_active)
        return;

    vec2 const &p1 = inputs->mouse.position.current;
    vec2 const &p0 = inputs->mouse.position.previous;

    bool const event_valid = !inputs->mouse.on_gui;
    bool const click_left = inputs->mouse.click.left;
    bool const ctrl = inputs->keyboard.ctrl;
    bool const shift = inputs->keyboard.shift;

    if (event_valid)
    { // If the mouse cursor is not on the ImGui area
        if (click_left && !ctrl && !shift) {
            camera_model.manipulator_translate_in_plane(p1 - p0);
        }
        if (click_left && ctrl && !shift) {
            if (dim == DIM_2D) {
                cam_projection.zoom(1.0f + (p1 - p0).y * 0.1f);
            }
            else {
                camera_model.manipulator_scale_distance_to_center((p1 - p0).y);
            }
        }
        if (click_left && shift && !ctrl) {
            if (dim == DIM_3D) {
                camera_model.manipulator_rotate_arcball(p0, p1);
            }
        }
    }

    update(camera_matrix_view);
}

void camera_controller::idle_frame(mat4 &camera_matrix_view)
{
    // Preconditions
    assert_cgp_no_msg(inputs != nullptr);
    assert_cgp_no_msg(window != nullptr);
    if (!is_active)
        return;

    float const angle_magnitude = inputs->time_interval;
    if (inputs->keyboard.left || inputs->keyboard.is_pressed(GLFW_KEY_R))
    {
        camera_model.manipulator_rotate_roll_pitch_yaw(angle_magnitude, 0, 0);
    }
    if (inputs->keyboard.right || inputs->keyboard.is_pressed(GLFW_KEY_F))
    {
        camera_model.manipulator_rotate_roll_pitch_yaw(-angle_magnitude, 0, 0);
    }

    update(camera_matrix_view);
}

void camera_controller::look_at(vec3 const &eye, vec3 const &center, vec3 const &up)
{
    camera_model.look_at(eye, center, up);
}

std::string camera_controller::doc_usage() const
{
    std::string doc;
    doc += "Info Camera Controller: Orbit - Camera that rotates around a central focus point.\n";
    doc += "Camera control: \n";
    doc += "   - Mouse left click + drag: ArcBall rotation.\n";
    doc += "   - Mouse right click + drag: Camera move close/far from its central focus point (the focus point remains unchanged).\n";
    doc += "   - Ctrl + Mouse left click + drag: Translate/Pan the camera and its central focus point in the viewspace plane.\n";
    doc += "   - Ctrl + Mouse right click + drag: Translate the camera and its central focus point in front/back direction.\n";
    doc += "   - Key left/right: Twist/Roll the camera angle around the view direction.\n";

    return doc;
}