#ifndef CAMERA_H
#define CAMERA_H

#include <glm/glm.hpp>

class Camera {
public:
	glm::mat4 get_view_matrix() const;
	glm::mat4 get_projection_matrix() const;
	glm::mat4 computeLookAt(glm::vec3 eye, glm::vec3 center, glm::vec3 up);
	// FIXME: add functions to manipulate camera objects.
	void rotate(glm::vec2 direction);
	void zoom(float direction);
	void forward_back(float direction);
	void strafe(float direction);
	void roll(float direction);
	void up_down(float direction);
	float mouse_x_ = 0.0f;
	float mouse_y_ = 0.0f;
private:
	//float camera_distance_ = 3.0;
	//glm::vec3 look_ = glm::vec3(0.0f, 0.0f, -1.0f);
	//glm::vec3 up_ = glm::vec3(0.0f, 1.0, 0.0f);
	//glm::vec3 eye_ = glm::vec3(0.0f, 0.0f, camera_distance_);
	glm::vec3 eye_ = glm::vec3(0.0f, 10.0f, 10.0f);
	glm::vec3 look_ = glm::normalize(glm::vec3(0.0f, 0.0f, 0.0f) - eye_);
	glm::vec3 up_ = glm::normalize(glm::cross(glm::vec3(1.0f,0.0f,0.0f), look_));
	float camera_distance_ = glm::length(eye_);
	// Note: you may need additional member variables
	glm::mat4 view_matrix = glm::mat4(1.0);
	glm::mat4 projection_matrix;
};

#endif
