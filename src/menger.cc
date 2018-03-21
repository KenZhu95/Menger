#include "menger.h"

namespace {
	const int kMinLevel = 0;
	const int kMaxLevel = 4;
};

Menger::Menger()
{
	// Add additional initialization if you like
}

Menger::~Menger()
{
}

void
Menger::set_nesting_level(int level)
{
	nesting_level_ = level;
	dirty_ = true;
}

bool
Menger::is_dirty() const
{
	return dirty_;
}

void
Menger::set_clean()
{
	dirty_ = false;
}

// FIXME generate Menger sponge geometry

void Menger::generate_geometry(std::vector<glm::vec4>& obj_vertices,
							   std::vector<glm::uvec3>& obj_faces) const
{
	if (!this->nesting_level_){
		create_menger(obj_vertices, obj_faces, min_bound, max_bound);
		return;
	}

	std::vector<glm::vec3> min_vertices;
	min_vertices.emplace_back(this->min_bound);

	for (int i = 0; i < nesting_level_; i++){
		glm::vec3 diff = (max_bound - min_bound) * float(1.0f / pow(3.0f, i+1));
		glm::vec3 x_diff = glm::vec3(diff.x, 0.0, 0.0);
		glm::vec3 y_diff = glm::vec3(0.0, diff.y, 0.0);
		glm::vec3 z_diff = glm::vec3(0.0, 0.0, diff.z);

		std::vector<glm::vec3> sub_cubes_min;
		sub_cubes_min = min_vertices;
		min_vertices.clear();

		for (glm::vec3 sub_cube : sub_cubes_min){
			min_vertices.emplace_back(sub_cube);
			min_vertices.emplace_back(sub_cube + x_diff);
			min_vertices.emplace_back(sub_cube + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + z_diff);
			min_vertices.emplace_back(sub_cube + z_diff + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + z_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + z_diff * 2.0f + x_diff);
			min_vertices.emplace_back(sub_cube + z_diff * 2.0f + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff);
			min_vertices.emplace_back(sub_cube + y_diff + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff + z_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff + z_diff * 2.0f + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + x_diff);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + z_diff);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + z_diff + x_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + z_diff * 2.0f);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + z_diff * 2.0f + x_diff);
			min_vertices.emplace_back(sub_cube + y_diff * 2.0f + z_diff * 2.0f + x_diff * 2.0f);

		}
	}

	glm::vec3 diff_min = (max_bound - min_bound) * float(1.0f / pow(3.0f, nesting_level_));
	for (auto it = min_vertices.begin(); it != min_vertices.end(); it++){
		create_menger(obj_vertices, obj_faces, *it, *it + diff_min);
	}

}

void
Menger::create_menger(std::vector<glm::vec4>& obj_vertices,
					  std::vector<glm::uvec3>& obj_faces,
					  glm::vec3 min,
					  glm::vec3 max) const
{
	unsigned long v_num = obj_vertices.size();

	obj_vertices.emplace_back(glm::vec4(min[0], min[1], min[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(max[0], min[1], min[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(min[0], max[1], min[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(max[0], max[1], min[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(max[0], min[1], max[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(min[0], min[1], max[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(max[0], max[1], max[2], 1.0f));
	obj_vertices.emplace_back(glm::vec4(min[0], max[1], max[2], 1.0f));


	obj_faces.emplace_back(glm::uvec3(v_num, v_num + 1, v_num + 2));
    obj_faces.emplace_back(glm::uvec3(v_num + 3, v_num + 2, v_num + 1));

	obj_faces.emplace_back(glm::uvec3(v_num + 1, v_num + 4, v_num + 3));
	obj_faces.emplace_back(glm::uvec3(v_num + 6, v_num + 3, v_num + 4));

	obj_faces.emplace_back(glm::uvec3(v_num + 5, v_num + 4, v_num + 7));
	obj_faces.emplace_back(glm::uvec3(v_num + 6, v_num + 7, v_num + 4));

	obj_faces.emplace_back(glm::uvec3(v_num + 0, v_num + 5, v_num + 2));
	obj_faces.emplace_back(glm::uvec3(v_num + 7, v_num + 2, v_num + 5));

	obj_faces.emplace_back(glm::uvec3(v_num + 0, v_num + 1, v_num + 5));
	obj_faces.emplace_back(glm::uvec3(v_num + 4, v_num + 5, v_num + 1));

	obj_faces.emplace_back(glm::uvec3(v_num + 2, v_num + 3, v_num + 7));
	obj_faces.emplace_back(glm::uvec3(v_num + 6, v_num + 7, v_num + 3));
}