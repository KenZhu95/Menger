#include "camera.h"
#include <glm/gtx/rotate_vector.hpp>

using namespace std;

namespace {
	float pan_speed = 0.1f;
	float roll_speed = 0.1f;
	float rotation_speed = 0.05f;
	float zoom_speed = 0.1f;
};

// FIXME: Calculate the view matrix
glm::mat4 Camera::get_view_matrix() const
{

	glm::vec3 center = eye_ + camera_distance_ * look_;
	glm::vec3 eye = eye_;
	glm::vec3 up = up_;

	glm::vec3 look = glm::normalize(center - eye);
	glm::vec3 tangent = glm::normalize(glm::cross(look, up));
	glm::vec3 up_n = glm::normalize(up);

	glm::mat4 result = glm::mat4(1.0);
	result[0][0] = tangent.x;
	result[1][0] = tangent.y;
	result[2][0] = tangent.z;
	result[0][1] = up_n.x;
	result[1][1] = up_n.y;
	result[2][1] = up_n.z;
	result[0][2] = -look.x;
	result[1][2] = -look.y;
	result[2][2] = -look.z;
	result[3][0] = -glm::dot(tangent, eye);
	result[3][1] = -glm::dot(up_n, eye);
	result[3][2] = glm::dot(look, eye);
	return result;
}

glm::mat4 Camera::get_projection_matrix() const{
	return projection_matrix;
}

glm::mat4 Camera::computeLookAt(glm::vec3 eye, glm::vec3 center, glm::vec3 up){
	glm::vec3 look = glm::normalize(center - eye);
	glm::vec3 tangent = glm::normalize(glm::cross(look, up));
	glm::vec3 up_n = glm::normalize(up);

	glm::mat4 result = glm::mat4(1.0);
	result[0][0] = tangent.x;
	result[1][0] = tangent.y;
	result[2][0] = tangent.z;
	result[0][1] = up_n.x;
	result[1][1] = up_n.y;
	result[2][1] = up_n.z;
	result[0][2] = -look.x;
	result[1][2] = -look.y;
	result[2][2] = -look.z;
	result[3][0] = -glm::dot(tangent, eye);
	result[3][1] = -glm::dot(up_n, eye);
	result[3][2] = glm::dot(look, eye);

	return result;
}


void Camera::rotate(glm::vec2 direction){
    glm::vec3 center = eye_ + camera_distance_ * look_;
    glm::vec3 tangent = glm::normalize(glm::cross(look_, up_));
    glm::vec3 world = up_ * direction.y + tangent * direction.x;
    glm::vec3 norm = glm::normalize(glm::cross(world, look_));
    look_ = glm::normalize(glm::rotate(look_, rotation_speed, norm));
    up_ = glm::normalize(glm::rotate(up_, -rotation_speed, norm));
}

void Camera::zoom(float direction){
	glm::vec3 center = eye_ + camera_distance_ * look_;
	camera_distance_ += zoom_speed * direction;
	eye_ = center - camera_distance_ * look_;
}

void Camera::forward_back(float direction){
	eye_ += direction * zoom_speed * glm::normalize(look_);
}

void Camera::strafe(float direction){
	glm::vec3 tangent = glm::normalize(glm::cross(look_, up_));
	eye_ += direction * pan_speed * tangent;
}

void Camera::roll(float direction){
	up_ = glm::rotate(up_, -roll_speed * direction, -look_);
}

void Camera::up_down(float direction){
	eye_ += direction * pan_speed * up_;
}


