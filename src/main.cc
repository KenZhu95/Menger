#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <ios>
#include <stdint.h>
#include <zlib.h>
#include <png.h>
#include "pngimage.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

// OpenGL library includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <debuggl.h>
#include "menger.h"
#include "camera.h"
#include <stdint.h>

#define FRONT "../../src/cubemap/negz.png"
#define BACK "../../src/cubemap/posz.png"
#define TOP "../../src/cubemap/posy.png"
#define BOTTOM "../../src/cubemap/negy.png"
#define LEFT "../../src/cubemap/negx.png"
#define RIGHT "../../src/cubemap/posx.png"

using std::string;



/* current versions of libpng should provide this macro: */
#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr)   ((png_ptr)->jmpbuf)
#endif

#ifdef DEBUG
#  define Trace(x)  {fprintf x ; fflush(stderr); fflush(stdout);}
#else
#  define Trace(x)  ;
#endif

void png_version_info(void)
{
	fprintf(stderr, "   Compiled with libpng %s; using libpng %s.\n",
		PNG_LIBPNG_VER_STRING, png_libpng_ver);
	fprintf(stderr, "   Compiled with zlib %s; using zlib %s.\n",
		ZLIB_VERSION, zlib_version);
}




namespace {

struct PNGReader {
	typedef unsigned char   uch;
	typedef unsigned short  ush;
	typedef unsigned long   ulg;

	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;

	png_uint_32  width, height;
	int  bit_depth, color_type;

	enum {
		SUCCESS = 0,
		BAD_SIGNATURE = 1,
		BAD_IHDR = 2,
		NO_MEMORY = 4,
		FOPEN_FAILURE = 8
	};

	~PNGReader()
	{
		if (png_ptr && info_ptr) {
			png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
			png_ptr = NULL;
			info_ptr = NULL;
		}
	}

	/*
 	 * return value:
	 * 0 for success,
	 * 1 for bad sig,
	 * 2 for bad IHDR,
	 * 4 for no mem,
	 * 8 for file open failure
	 */
	int init(const char* filename)
	{
		uch sig[8];
		FILE *infile;

		if ((infile = fopen(filename, "rb")) == NULL){
			return FOPEN_FAILURE;
        }
		/* check that the file really is a PNG image; could
		 * have used slightly more general png_sig_cmp() function instead */

		fread(sig, 1, 8, infile);
		if (png_sig_cmp(sig, 0, 8) != 0)
			return BAD_SIGNATURE;
		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if (!png_ptr)
			return NO_MEMORY;

		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr) {
			png_destroy_read_struct(&png_ptr, NULL, NULL);
			return NO_MEMORY;
		}

		/* we could create a second info struct here (end_info), but it's only
		 * useful if we want to keep pre- and post-IDAT chunk info separated
		 * (mainly for PNG-aware image editors and converters) */

		/* setjmp() must be called in every function that calls a PNG-reading
		 * libpng function */

		if (setjmp(png_jmpbuf(png_ptr))) {
			png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
			return BAD_IHDR;
		}

		png_init_io(png_ptr, infile);
		png_set_sig_bytes(png_ptr, 8);  /* we already read the 8 signature bytes */
		png_read_info(png_ptr, info_ptr);  /* read all PNG info up to image data */


		/* alternatively, could make separate calls to png_get_image_width(),
		 * etc., but want bit_depth and color_type for later [don't care about
		 * compression_type and filter_type => NULLs] */

		png_get_IHDR(png_ptr, info_ptr,
				&width, &height, &bit_depth, &color_type,
				NULL, NULL, NULL);
		/* OK, that's all we need for now; return happy */

		return 0;
	}

	uch *get_image(double display_exponent, int &pChannels, int &pRowbytes)
	{
		double  gamma;
		png_uint_32  i, rowbytes;

		/* setjmp() must be called in every function that calls a PNG-reading
		 * libpng function */

		if (setjmp(png_jmpbuf(png_ptr))) {
			png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
			return NULL;
		}

		/* expand palette images to RGB, low-bit-depth grayscale images to 8 bits,
		 * transparency chunks to full alpha channel; strip 16-bit-per-sample
		 * images to 8 bits per sample; and convert grayscale to RGB[A] */

		if (color_type == PNG_COLOR_TYPE_PALETTE)
			png_set_expand(png_ptr);
		if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
			png_set_expand(png_ptr);
		if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
			png_set_expand(png_ptr);
		if (bit_depth == 16)
			png_set_strip_16(png_ptr);
		if (color_type == PNG_COLOR_TYPE_GRAY ||
				color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
			png_set_gray_to_rgb(png_ptr);

		/* unlike the example in the libpng documentation, we have *no* idea where
		 * this file may have come from--so if it doesn't have a file gamma, don't
		 * do any correction ("do no harm") */

		if (png_get_gAMA(png_ptr, info_ptr, &gamma))
			png_set_gamma(png_ptr, display_exponent, gamma);

		/* all transformations have been registered; now update info_ptr data,
		 * get rowbytes and channels, and allocate image memory */

		png_read_update_info(png_ptr, info_ptr);

		pRowbytes = rowbytes = png_get_rowbytes(png_ptr, info_ptr);
		pChannels = (int)png_get_channels(png_ptr, info_ptr);

		uch *image_data = NULL;

		if ((image_data = (uch *)malloc(rowbytes*height)) == NULL) {
			png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
			return NULL;
		}

		std::vector<png_bytep> row_pointers(height);

		Trace((stderr, "readpng_get_image:  channels = %d, rowbytes = %ld, height = %ld\n", *pChannels, rowbytes, height));

		/* set the individual row_pointers to point at the correct offsets */

		for (i = 0;  i < height;  ++i)
			row_pointers[i] = image_data + i*rowbytes;

		/* now we can go ahead and just read the whole image */

		png_read_image(png_ptr, row_pointers.data());

		/* and we're done!  (png_read_end() can be omitted if no processing of
		 * post-IDAT text/time/etc. is desired) */

		png_read_end(png_ptr, NULL);
		return image_data;
	}
};

}; // Anonymous namespace

unsigned char* readPNGData(const char *fname, int& width, int& height);

unsigned char* readPNGData(const char *fname, int& width, int& height){
	PNGReader reader;
	unsigned char* data = nullptr;
	if (reader.init(fname) != 0){
		return data; 
	}

	width = reader.width;
	height = reader.height;

	static constexpr double gamma = 2.2;
	int channels, rowBytes;
	unsigned char* indata = reader.get_image(gamma, channels, rowBytes);

	return indata;
}


int window_width = 800, window_height = 600;

// VBO and VAO descriptors.
enum { kVertexBuffer, kIndexBuffer, kNumVbos };

// These are our VAOs.
enum { kGeometryVao, kFloorVao, kOceanVao, kSphereVao, kNumVaos };

GLuint g_array_objects[kNumVaos];  // This will store the VAO descriptors.
GLuint g_buffer_objects[kNumVaos][kNumVbos];  // These will store VBO descriptors.


/* big cube. returns Vertex Array Object */
GLuint make_big_cube() {
	float points[] = {
		-20.0f, 20.0f,	-20.0f, -20.0f, -20.0f, -20.0f, 20.0f,	-20.0f, -20.0f,
		20.0f,	-20.0f, -20.0f, 20.0f,	20.0f,	-20.0f, -20.0f, 20.0f,	-20.0f,

		-20.0f, -20.0f, 20.0f,	-20.0f, -20.0f, -20.0f, -20.0f, 20.0f,	-20.0f,
		-20.0f, 20.0f,	-20.0f, -20.0f, 20.0f,	20.0f,	-20.0f, -20.0f, 20.0f,

		20.0f,	-20.0f, -20.0f, 20.0f,	-20.0f, 20.0f,	20.0f,	20.0f,	20.0f,
		20.0f,	20.0f,	20.0f,	20.0f,	20.0f,	-20.0f, 20.0f,	-20.0f, -20.0f,

		-20.0f, -20.0f, 20.0f,	-20.0f, 20.0f,	20.0f,	20.0f,	20.0f,	20.0f,
		20.0f,	20.0f,	20.0f,	20.0f,	-20.0f, 20.0f,	-20.0f, -20.0f, 20.0f,

		-20.0f, 20.0f,	-20.0f, 20.0f,	20.0f,	-20.0f, 20.0f,	20.0f,	20.0f,
		20.0f,	20.0f,	20.0f,	-20.0f, 20.0f,	20.0f,	-20.0f, 20.0f,	-20.0f,

		-20.0f, -20.0f, -20.0f, -20.0f, -20.0f, 20.0f,	20.0f,	-20.0f, -20.0f,
		20.0f,	-20.0f, -20.0f, -20.0f, -20.0f, 20.0f,	20.0f,	-20.0f, 20.0f
	};
	GLuint vbo;
	glGenBuffers( 1, &vbo );
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	glBufferData( GL_ARRAY_BUFFER, 3 * 36 * sizeof( GLfloat ), &points,
								GL_STATIC_DRAW );

	GLuint vao;
	glGenVertexArrays( 1, &vao );
	glBindVertexArray( vao );
	glEnableVertexAttribArray( 0 );
	glBindBuffer( GL_ARRAY_BUFFER, vbo );
	glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, NULL );
	return vao;
}

/* use stb_image to load an image file into memory, and then into one side of
a cube-map texture. */
bool load_cube_map_side( GLuint texture, GLenum side_target,
												 const char *file_name ) {
	glBindTexture( GL_TEXTURE_CUBE_MAP, texture );

	int x, y, n;
	int force_channels = 4;
	//unsigned char *image_data = stbi_load( file_name, &x, &y, &n, force_channels );
	unsigned char *image_data = readPNGData(file_name, x, y);

	if ( !image_data ) {
		fprintf( stderr, "ERROR: could not load %s\n", file_name );
		return false;
	}
	// non-power-of-2 dimensions check
	if ( ( x & ( x - 1 ) ) != 0 || ( y & ( y - 1 ) ) != 0 ) {
		fprintf( stderr, "WARNING: image %s is not power-of-2 dimensions\n",
						 file_name );
	}

	// copy image data into 'target' side of cube map
	glTexImage2D( side_target, 0, GL_RGB, x, y, 0, GL_RGB, GL_UNSIGNED_BYTE,
								image_data );
	free( image_data );
	return true;
}

/* load all 6 sides of the cube-map from images, then apply formatting to the
final texture */
void create_cube_map( const char *front, const char *back, const char *top,
											const char *bottom, const char *left, const char *right,
											GLuint *tex_cube ) {
	// generate a cube-map texture to hold all the sides
	glActiveTexture( GL_TEXTURE0 );
	glGenTextures( 1, tex_cube );

	// load each image and copy into a side of the cube-map texture
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, front ) );
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_POSITIVE_Z, back ) );
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, top ) );
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, bottom ) );
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, left ) );
	( load_cube_map_side( *tex_cube, GL_TEXTURE_CUBE_MAP_POSITIVE_X, right ) );
	// format cube map texture
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
}

// C++ 11 String Literal
// See http://en.cppreference.com/w/cpp/language/string_literal
const char* vertex_shader =
R"zzz(#version 400 core
in vec4 vertex_position;
uniform mat4 view;
uniform vec4 light_position;
out vec4 vs_light_direction;
out vec3 v_position;
void main()
{
	//gl_Position = view * vertex_position;
	//vs_light_direction = -gl_Position + view * light_position;
	gl_Position = vertex_position;
	vs_light_direction = light_position - view * vertex_position;
	v_position = vertex_position.xyz;
}
)zzz";


const char* geometry_shader =
R"zzz(#version 400 core
layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;
uniform mat4 projection;
uniform mat4 view;
in vec4 vs_light_direction[];
flat out vec4 normal;
out vec4 light_direction;
out vec4 world_position;
out vec3 bary_coor;
flat out vec3 heights;
void main()
{
	int n = 0;
	//normal = vec4(0.0, 0.0, 1.0f, 0.0);
	normal = vec4(normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz)),1.0);
	vec3 pos1 = gl_in[0].gl_Position.xyz;
	vec3 pos2 = gl_in[1].gl_Position.xyz;
	vec3 pos3 = gl_in[2].gl_Position.xyz;
	vec3 l_a = pos2 - pos1;
	vec3 l_b = pos3 - pos1;
	float area = length(cross(pos2-pos1,pos3-pos1));
	float length1 = area / length(pos3-pos2);
	float length2 = area / length(pos3-pos1);
	float length3 = area / length(pos2-pos1);
	heights = vec3(length1, length2, length3);

	for (n = 0; n < gl_in.length(); n++) {
		light_direction = vs_light_direction[n];
		world_position = gl_in[n].gl_Position;
		//gl_Position = projection * gl_in[n].gl_Position;
		gl_Position = projection * view * gl_in[n].gl_Position;
		if (n == 0){
			bary_coor = vec3(1.0,0.0,0.0);
		} else if (n == 1){
			bary_coor = vec3(0.0, 1.0, 0.0);
		} else {
			bary_coor = vec3(0.0, 0.0, 1.0);
		}
		EmitVertex();
	}
	EndPrimitive();
}
)zzz";

const char* fragment_shader =
R"zzz(#version 400 core
flat in vec4 normal;
in vec4 light_direction;
out vec4 fragment_color;

struct normalColor{
	vec3 normal;
	vec4 color;
};

void main()
{
	vec4 color = vec4(0.5, 0.5, 0.5, 1.0);
	color = vec4(abs(normal.x), abs(normal.y), abs(normal.z), 0.0);
	float dot_nl = dot(normalize(light_direction), normalize(normal));
	dot_nl = clamp(dot_nl, 0.0, 1.0);
	//fragment_color = clamp(dot_nl * color, 0.0, 1.0);
	fragment_color = color;

}
)zzz";

const char* sphere_fragment_shader =
R"zzz(#version 400 core
flat in vec4 normal;
in vec4 light_direction;
out vec4 fragment_color;

struct normalColor{
	vec3 normal;
	vec4 color;
};

void main()
{
	fragment_color = vec4(1.0,1.0,1.0,0.0);

}
)zzz";


// FIXME: Implement shader effects with an alternative shader.
const char* floor_fragment_shader =
R"zzz(#version 400 core
in vec4 normal;
in vec4 light_direction;
in vec4 world_position;
in vec3 bary_coor;
flat in vec3 heights;
uniform vec4 light_position;
uniform bool if_render_frame;
out vec4 fragment_color;

void main()
{
	vec4 pos = world_position;
	float check_width = 1.0;
	float i = floor(pos.x / check_width);
	float j = floor(pos.z / check_width);
	vec3 color = mod(i+j, 2) * vec3(1.0, 1.0, 1.0);
	float epsilon = 0.02;
	if (if_render_frame && ((bary_coor.x * heights.x) < epsilon  || (bary_coor.y * heights.y) < epsilon  || (bary_coor.z * heights.z) < epsilon )){
		color = vec3(0.0, 1.0, 0.0);
	}
	fragment_color = vec4(color, 1.0);
	float dot_nl = dot(normalize(light_position - pos/pos.w), vec4(0,1.0f,0,1.0f));
	dot_nl = clamp(dot_nl,0.0,1.0);
	fragment_color = clamp(dot_nl * fragment_color, 0.0, 1.0);
	//fragment_color = vec4(1.0,1.0,0.0,1.0);

}
)zzz";

const char* floor_tess_control_shader =
R"zzz(#version 400 core
layout (vertices = 3) out;
in vec4 vs_light_direction[];
uniform float tessellation_level_inner;
uniform float tessellation_level_outer;
out vec4 tess_light_direction[];

void main(){
	gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
	tess_light_direction[gl_InvocationID] = vs_light_direction[gl_InvocationID];
	if (gl_InvocationID == 0){
		gl_TessLevelInner[0] = tessellation_level_inner;
		gl_TessLevelOuter[0] = tessellation_level_outer;
		gl_TessLevelOuter[1] = tessellation_level_outer;
		gl_TessLevelOuter[2] = tessellation_level_outer;
	}
}
)zzz";

const char* floor_tess_evaluation_shader =
R"zzz(#version 400 core
layout(triangles) in;
in vec4 tess_light_direction[];
out vec4 vs_light_direction;

void main(){
	gl_Position = (gl_TessCoord.x * gl_in[0].gl_Position) + (gl_TessCoord.y * gl_in[1].gl_Position) + (gl_TessCoord.z * gl_in[2].gl_Position);
	vs_light_direction = (gl_TessCoord.x * tess_light_direction[0]) + (gl_TessCoord.y * tess_light_direction[1]) + (gl_TessCoord.z * tess_light_direction[2]);
}
)zzz";


const char* ocean_tess_control_shader =
R"zzz(#version 400 core
layout (vertices = 4) out;
in vec4 vs_light_direction[];
in vec3 v_position[];
uniform float tidal_time;
uniform bool if_generate_tidal; 
uniform float tessellation_level_inner;
uniform float tessellation_level_outer;
out vec4 tess_light_direction[];

void main(){
	gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
	tess_light_direction[gl_InvocationID] = vs_light_direction[gl_InvocationID];

    vec2 tidal_centre;
	if (if_generate_tidal){
    	tidal_centre = vec2(0.5 * tidal_time, 0.0);
    } else {
    	tidal_centre = vec2(0.0, 0.0);
    }

    vec3 center = vec3(tidal_centre.x, -2.0, tidal_centre.y);
    float min_distance = min(min(length(v_position[0] - center), length(v_position[1] - center)), min(length(v_position[2] - center), length(v_position[3] - center)));
    float multi = clamp(4 - min_distance, 1.0, 4.0);

	if (gl_InvocationID == 0){
		gl_TessLevelInner[0] = tessellation_level_inner;
		gl_TessLevelInner[1] = tessellation_level_inner;
		gl_TessLevelOuter[0] = tessellation_level_outer;
		gl_TessLevelOuter[1] = tessellation_level_outer;
		gl_TessLevelOuter[2] = tessellation_level_outer;
		gl_TessLevelOuter[3] = tessellation_level_outer;
	}

	if (if_generate_tidal) {
		gl_TessLevelInner[0] = tessellation_level_inner * multi;
		gl_TessLevelInner[1] = tessellation_level_inner * multi;
		gl_TessLevelOuter[0] = tessellation_level_outer * multi;
		gl_TessLevelOuter[1] = tessellation_level_outer * multi;
		gl_TessLevelOuter[2] = tessellation_level_outer * multi;
		gl_TessLevelOuter[3] = tessellation_level_outer * multi;
	}


}
)zzz";


const char* ocean_tess_evaluation_shader =
R"zzz(#version 400 core
layout(quads) in;
in vec4 tess_light_direction[];
out vec4 vs_light_direction;

void main(){
	vec4 cor_x = mix(gl_in[0].gl_Position, gl_in[1].gl_Position, gl_TessCoord.x);
	vec4 cor_y = mix(gl_in[3].gl_Position, gl_in[2].gl_Position, gl_TessCoord.x);
	gl_Position = mix(cor_x, cor_y, gl_TessCoord.y);
}
)zzz";


const char* ocean_fragment_shader = 
R"zzz(#version 400 core
in vec4 normal;
in vec4 light_direction;
in vec4 world_position;
in vec3 bary_coor;
flat in vec4[3] new_position;
flat in vec4[3] old_position;
uniform vec4 light_position;
uniform bool if_render_frame;
uniform mat4 view;
out vec4 fragment_color;

void main()
{
	vec3 pos1 = new_position[0].xyz;
	vec3 pos2 = new_position[1].xyz;
	vec3 pos3 = new_position[2].xyz;
	float area = length(cross(pos2-pos1,pos3-pos1));
	float length1 = area / length(pos3-pos2);
	float length2 = area / length(pos3-pos1);
	float length3 = area / length(pos2-pos1);
	vec3 heights = vec3(length1, length2, length3);
	vec3 barys;
	float area_1 = length(cross(world_position.xyz-pos2, world_position.xyz-pos3));
	float area_2 = length(cross(world_position.xyz-pos1, world_position.xyz-pos3));
	float area_3 = length(cross(world_position.xyz-pos2, world_position.xyz-pos1));
	barys.x = area_1 / area;
	barys.y = area_2 / area;
	barys.z = area_3 / area;

	vec4 pos = world_position;
	vec3 color = vec3(0.0, 0.4, 0.6);
	float epsilon = 0.02;
	if (if_render_frame && ((barys.x * heights.x) < epsilon  || (barys.y * heights.y) < epsilon  || (barys.z * heights.z) < epsilon )){
		color = vec3(0.0, 1.0, 0.0);
		fragment_color = vec4(color, 1.0);
	} else {
		fragment_color = vec4(color, 1.0);
		//diffuse term
		float dot_nl = dot(normalize(light_position - pos/pos.w), normalize(normal));
		//float dot_nl = dot(normalize(light_direction), normalize(normal));
		dot_nl = clamp(dot_nl,0.0,1.0);
		fragment_color = clamp(dot_nl * fragment_color, 0.0, 1.0);
		//specular term
		vec4 vector_r = normalize(reflect(pos / pos.w -light_position, normal));
		//vec4 vector_r = normalize(reflect(-light_direction, normal));
		vec4 vector_v = normalize((inverse(view)[3] - pos));
		vec4 specular = vec4(color, 1.0) * pow(clamp(max(dot(vector_r, vector_v),0),0.0,1.0), 3);
		fragment_color += specular;
		fragment_color = clamp(fragment_color, 0.0, 1.0);

	}


}

)zzz";


const char* ocean_geometry_shader = 
R"zzz(#version 400 core
layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;
uniform mat4 projection;
uniform mat4 view;
uniform float time;
uniform float tidal_time;
uniform bool if_generate_tidal;
in vec4 vs_light_direction[];
out vec4 normal;
out vec4 light_direction;
out vec4 world_position;
out vec3 bary_coor;
flat out vec4[3] new_position;
flat out vec4[3] old_position;
flat out vec2 tidal_center;

float wave_height(float amplitude, vec2 direction, vec2 coor_xz, float frequency, float time, float speed){
	return amplitude * sin(dot(direction, coor_xz) * frequency + time * speed * frequency);
}

vec4 wave_normal(float amplitude, vec2 direction, vec2 coor_xz, float frequency, float time, float speed){
	float normal_x = amplitude * frequency * cos(dot(direction, coor_xz) * frequency + time * speed * frequency) * direction.x;
	float normal_z = amplitude * frequency * cos(dot(direction, coor_xz) * frequency + time * speed * frequency) * direction.y;
	return vec4(normal_x, 1.0f, normal_z, 0.0f);
}

void main()
{
	int n = 0;
	//normal = vec4(0.0, 0.0, 1.0f, 0.0);
	//normal = vec4(normalize(cross(gl_in[1].gl_Position.xyz - gl_in[0].gl_Position.xyz, gl_in[2].gl_Position.xyz - gl_in[0].gl_Position.xyz)),1.0);



	//rquire to recalculate y coordinate for world_position

    //calculate Gaussian center

    if (if_generate_tidal){
    	tidal_center = vec2(0.5 * tidal_time, 0.0);
    } else {
    	tidal_center = vec2(0.0, 0.0);
    }


	for (n = 0; n < gl_in.length(); n++) {
		light_direction = vs_light_direction[n];
		vec4 temp_world = gl_in[n].gl_Position;
		old_position[n] = gl_in[n].gl_Position;
		temp_world.y += wave_height(0.4f, vec2(1.0f, 0.0f), vec2(temp_world.x, temp_world.z), 1.0f, time, 1.0f);
		temp_world.y += wave_height(0.8f, vec2(0.0f, 1.0f), vec2(temp_world.x, temp_world.z), 0.5f, time, 0.8f);
		temp_world.y += wave_height(0.5f, normalize(vec2(1.0f, 1.0f)), vec2(temp_world.x, temp_world.z), 0.6f, time, 2.0f);
		temp_world.y += wave_height(0.3f, normalize(vec2(-1.0f, 1.0f)), vec2(temp_world.x, temp_world.z), 2.0f, time, 1.0f);
		if (if_generate_tidal){
			temp_world.y += 10 * exp(-2*pow((temp_world.x - tidal_center.x),2)-2*pow((temp_world.z),2));
		}
		world_position =temp_world;
		//world_position = gl_in[n].gl_Position;
		new_position[n] = temp_world;
		normal = wave_normal(0.4f, vec2(1.0f, 0.0f), vec2(temp_world.x, temp_world.z), 1.0f, time, 1.0f);
		normal += wave_normal(0.8f, vec2(0.0f, 1.0f), vec2(temp_world.x, temp_world.z), 0.5f, time, 0.8f);
		normal += wave_normal(0.5f, normalize(vec2(1.0f, 1.0f)), vec2(temp_world.x, temp_world.z), 0.6f, time, 2.0f);
		normal += wave_normal(0.3f, normalize(vec2(-1.0f, 1.0f)), vec2(temp_world.x, temp_world.z), 2.0f, time, 1.0f);
		normal = normalize(normal);
		//gl_Position = projection * gl_in[n].gl_Position;
		gl_Position = projection * view * temp_world;
		if (n == 0){
			bary_coor = vec3(1.0,0.0,0.0);
		} else if (n == 1){
			bary_coor = vec3(0.0, 1.0, 0.0);
		} else {
			bary_coor = vec3(0.0, 0.0, 1.0);
		}
		EmitVertex();
	}
	EndPrimitive();


}
)zzz";

const char* skybox_fragment_shader = 
R"zzz(#version 410 core
in vec3 texcoords;
uniform samplerCube cube_texture;
out vec4 frag_colour;

void main () {
	frag_colour = texture (cube_texture, texcoords);
}

)zzz";


const char* skybox_vertex_shader = 
R"zzz(#version 400 core
in vec3 vp;
uniform mat4 projection, view;
out vec3 texcoords;

void main () {
	texcoords = vp;
	gl_Position = projection * view * vec4 (vp, 1.0);
}
)zzz";



void
CreateTriangle(std::vector<glm::vec4>& vertices,
        std::vector<glm::uvec3>& indices)
{
	vertices.push_back(glm::vec4(-0.5f, -0.5f, -0.5f, 1.0f));
	vertices.push_back(glm::vec4(0.5f, -0.5f, -0.5f, 1.0f));
	vertices.push_back(glm::vec4(0.0f, 0.5f, -0.5f, 1.0f));
	indices.push_back(glm::uvec3(0, 1, 2));
}

// FIXME: Save geometry to OBJ file
void
SaveObj(const std::string& file,
        const std::vector<glm::vec4>& vertices,
        const std::vector<glm::uvec3>& indices)
{
	std::ofstream f;
	f.open(file);
	for (auto &vertex : vertices){
		f << "v" << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
	}

	for (auto &index : indices){
		f << "f" << (index.x + 1) << " " << (index.y + 1) << " " << (index.z + 1) << std::endl;
	}

	f.close();
}

void
ErrorCallback(int error, const char* description)
{
	std::cerr << "GLFW Error: " << description << "\n";
}

std::shared_ptr<Menger> g_menger;
Camera g_camera;
bool FPSMode = true;
bool if_save = false;
bool if_faces_render = true;
bool if_ocean_mode = false;
bool if_generate_tidal =false;
bool if_reset_tidal_time = false;
bool if_render_frame = true;
float tessellation_level_outer = 4.0f;
float tessellation_level_inner = 4.0f;

void
KeyCallback(GLFWwindow* window,
            int key,
            int scancode,
            int action,
            int mods)
{
	// Note:
	// This is only a list of functions to implement.
	// you may want to re-organize this piece of code.
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
	else if (key == GLFW_KEY_S && mods == GLFW_MOD_CONTROL && action == GLFW_RELEASE) {
		// FIXME: save geometry to OBJ
		if_save = true;
		//SaveObj("geometry.obj", obj_vertices, obj_faces);
		
	} else if (key == GLFW_KEY_F && mods == GLFW_MOD_CONTROL && action == GLFW_RELEASE) {
		if_faces_render = !if_faces_render;		
	} else if (key == GLFW_KEY_O && mods == GLFW_MOD_CONTROL && action == GLFW_RELEASE) {
		if_ocean_mode = !if_ocean_mode;		
	} else if (key == GLFW_KEY_T && mods == GLFW_MOD_CONTROL && action == GLFW_RELEASE) {
		if_generate_tidal = true;
		if_reset_tidal_time = false;		
	} else if (key == GLFW_KEY_F && action != GLFW_RELEASE && mods!= GLFW_MOD_CONTROL) {
		if_render_frame = !if_render_frame;
	} else if (key == GLFW_KEY_W && action != GLFW_RELEASE) {
		// FIXME: WASD
		g_camera.forward_back(1.0);
	} else if (key == GLFW_KEY_S && action != GLFW_RELEASE && mods!= GLFW_MOD_CONTROL) {
		g_camera.forward_back(-1.0);
	} else if (key == GLFW_KEY_A && action != GLFW_RELEASE) {
		g_camera.strafe(-1.0);
	} else if (key == GLFW_KEY_D && action != GLFW_RELEASE) {
		g_camera.strafe(1.0);
	} else if (key == GLFW_KEY_LEFT && action != GLFW_RELEASE) {
		// FIXME: Left Right Up and Down
		g_camera.roll(-1.0);
	} else if (key == GLFW_KEY_RIGHT && action != GLFW_RELEASE) {
		g_camera.roll(1.0);
	} else if (key == GLFW_KEY_DOWN && action != GLFW_RELEASE) {
		g_camera.up_down(-1.0);
	} else if (key == GLFW_KEY_UP && action != GLFW_RELEASE) {
		g_camera.up_down(1.0);
	} else if (key == GLFW_KEY_C && action != GLFW_RELEASE) {
		// FIXME: FPS mode on/off
		FPSMode = !FPSMode;
	} else if (key == GLFW_KEY_MINUS && action != GLFW_RELEASE) {
		if (tessellation_level_outer > 1.0f){
			tessellation_level_outer -= 1.0f;
		}
	} else if (key == GLFW_KEY_EQUAL && action != GLFW_RELEASE) {
		tessellation_level_outer += 1.0f;
	} else if (key == GLFW_KEY_COMMA && action != GLFW_RELEASE) {
		if (tessellation_level_inner > 1.0f) {
			tessellation_level_inner -= 1.0f;
		}
	} else if (key == GLFW_KEY_PERIOD && action != GLFW_RELEASE) {
		tessellation_level_inner += 1.0f;
	}
	if (!g_menger)
		return ; // 0-4 only available in Menger mode.
	if (key == GLFW_KEY_0 && action != GLFW_RELEASE) {
		// FIXME: Change nesting level of g_menger
		// Note: GLFW_KEY_0 - 4 may not be continuous.
		g_menger->set_nesting_level(0);
	} else if (key == GLFW_KEY_1 && action != GLFW_RELEASE) {
		g_menger->set_nesting_level(1);
	} else if (key == GLFW_KEY_2 && action != GLFW_RELEASE) {
		g_menger->set_nesting_level(2);
	} else if (key == GLFW_KEY_3 && action != GLFW_RELEASE) {
		g_menger->set_nesting_level(3);
	} else if (key == GLFW_KEY_4 && action != GLFW_RELEASE) {
		g_menger->set_nesting_level(4);
	}
}

int g_current_button;
bool g_mouse_pressed;
bool g_prev_pressed;

void generate_floor(std::vector<glm::vec4>& vertices, std::vector<glm::uvec3>& faces);
void generate_ocean(std::vector<glm::vec4>& vertices, std::vector<glm::uvec4>& faces);
void generate_sphere(std::vector<glm::vec4>& vertices, std::vector<glm::uvec3>& faces);

void
MousePosCallback(GLFWwindow* window, double mouse_x, double mouse_y)
{
	glm::vec2 mouse(mouse_x, mouse_y);
	glm::vec2 prev(g_camera.mouse_x_, g_camera.mouse_y_);
	glm::vec2 diff = mouse - prev;
	g_prev_pressed = g_mouse_pressed;
	g_camera.mouse_x_ = mouse_x;
	g_camera.mouse_y_ = mouse_y;
	if (!g_mouse_pressed)
		return;
	if (!g_prev_pressed){
		//std::cout << "not pressed previously" << std::endl;
		return;
	}
	if (g_current_button == GLFW_MOUSE_BUTTON_LEFT) {
		// FIXME: left drag
		g_camera.rotate(glm::vec2(diff.x, diff.y));
	} else if (g_current_button == GLFW_MOUSE_BUTTON_RIGHT) {
		// FIXME: right drag
		//std::cout << "Now: " << mouse_y << std::endl;
		g_camera.zoom(10.0f * diff.y / window_height);
	} else if (g_current_button == GLFW_MOUSE_BUTTON_MIDDLE) {
		// FIXME: middle drag

	}

	//g_camera.mouse_x_ = mouse_x;
	//g_camera.mouse_y_ = mouse_y;
	//g_prev_pressed = g_mouse_pressed;
}

void
MouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	g_mouse_pressed = (action == GLFW_PRESS);
	g_current_button = button;
}

void generate_floor(std::vector<glm::vec4>& vertices, std::vector<glm::uvec3>& faces){
	unsigned long v_num = vertices.size();
	vertices.emplace_back(glm::vec4(10, -2, 10, 1));
	vertices.emplace_back(glm::vec4(10, -2, -10, 1));
	vertices.emplace_back(glm::vec4(-10, -2, -10, 1));
	vertices.emplace_back(glm::vec4(-10, -2, 10, 1));

	faces.emplace_back(glm::uvec3(v_num, v_num + 1, v_num + 2));
	faces.emplace_back(glm::uvec3(v_num, v_num + 2, v_num + 3));
}

void generate_ocean(std::vector<glm::vec4>& vertices, std::vector<glm::uvec4>& faces){
	for (int i=0; i< 16; i++){
		for (int j=0; j<16; j++){
			unsigned long v_num = vertices.size();
			vertices.emplace_back(glm::vec4(-20.0f + 2.5f * i, -2.0f, -20.0f + 2.5f * j, 1.0f));
			vertices.emplace_back(glm::vec4(-20.0f + 2.5f * i, -2.0f, -20.0f + 2.5f * (j+1), 1.0f));
			vertices.emplace_back(glm::vec4(-20.0f + 2.5f * (i+1), -2.0f, -20.0f + 2.5f * j, 1.0f));
			vertices.emplace_back(glm::vec4(-20.0f + 2.5f * (i+1), -2.0f, -20.0f + 2.5f * (j+1), 1.0f));
			faces.emplace_back(glm::uvec4(v_num, v_num + 2, v_num + 3, v_num + 1));

		}
	}
}

const double pi = std::acos(-1);

void generate_sphere(std::vector<glm::vec4>& vertices, std::vector<glm::uvec3>& faces){
	int i_max = 720;
	int j_max = 360;
	float a = pi/360;
	for (int i=0; i<i_max; i++){
		for (int j=0; j<j_max; j++){
			unsigned long v_num = vertices.size();
			glm::vec4 light_position = glm::vec4(-10.0f, 10.0f, 0.0f, 1.0f);
			vertices.emplace_back(light_position + glm::normalize(glm::vec4(std::cos((j-90)*a)*std::cos(i*a), std::sin((j-90)*a), std::cos((j-90)*a)*std::sin(i*a), 0.0)));
			vertices.emplace_back(light_position + glm::normalize(glm::vec4(std::cos((j+1-90)*a)*std::cos(i*a),  std::sin((j+1-90)*a), std::cos((j+1-90)*a)*std::sin(i*a), 0.0)));
			vertices.emplace_back(light_position + glm::normalize(glm::vec4(std::cos((j-90)*a)*std::cos((i+1)*a),  std::sin((j-90)*a), std::cos((j-90)*a)*std::sin((i+1)*a), 0.0)));
			vertices.emplace_back(light_position + glm::normalize(glm::vec4(std::cos((j+1-90)*a)*std::cos((i+1)*a),  std::sin((j+1-90)*a),std::cos((j+1-90)*a)*std::sin((i+1)*a), 0.0)));
			faces.emplace_back(glm::uvec3(v_num + 0, v_num + 1, v_num + 2));
			//faces.emplace_back(glm::uvec3(v_num + 0, v_num + 1, v_num + 3));
			//faces.emplace_back(glm::uvec3(v_num + 0, v_num + 2, v_num + 3));
			faces.emplace_back(glm::uvec3(v_num + 3, v_num + 2, v_num + 1));
			//faces.emplace_back(glm::uvec4(v_num, v_num + 2, v_num + 3, v_num + 1));
			//std::cout << vertices[v_num + 3].y << std::endl;
		}
	}
}

int main(int argc, char* argv[])
{
	std::string window_title = "Menger";
	if (!glfwInit()) exit(EXIT_FAILURE);
	g_menger = std::make_shared<Menger>();
	glfwSetErrorCallback(ErrorCallback);

	// Ask an OpenGL 4.1 core profile context
	// It is required on OSX and non-NVIDIA Linux
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	GLFWwindow* window = glfwCreateWindow(window_width, window_height,
			&window_title[0], nullptr, nullptr);
	CHECK_SUCCESS(window != nullptr);
	glfwMakeContextCurrent(window);
	glewExperimental = GL_TRUE;

	CHECK_SUCCESS(glewInit() == GLEW_OK);
	glGetError();  // clear GLEW's error for it
	glfwSetKeyCallback(window, KeyCallback);
	glfwSetCursorPosCallback(window, MousePosCallback);
	glfwSetMouseButtonCallback(window, MouseButtonCallback);
	glfwSwapInterval(1);
	const GLubyte* renderer = glGetString(GL_RENDERER);  // get renderer string
	const GLubyte* version = glGetString(GL_VERSION);    // version as a string
	std::cout << "Renderer: " << renderer << "\n";
	std::cout << "OpenGL version supported:" << version << "\n";

	std::vector<glm::vec4> obj_vertices;
	std::vector<glm::uvec3> obj_faces;

    //FIXME: Create the geometry from a Menger object.
    //CreateTriangle(obj_vertices, obj_faces);

	g_menger->set_nesting_level(4);
	g_menger->generate_geometry(obj_vertices, obj_faces);
	//g_menger->set_clean();

	std::vector<glm::vec4> floor_v;
	std::vector<glm::uvec3> floor_f;
	std::vector<glm::vec4> ocean_v;
	std::vector<glm::uvec4> ocean_f;
	std::vector<glm::vec4> sphere_v;
	std::vector<glm::uvec3> sphere_f;

	generate_floor(floor_v, floor_f);
	generate_ocean(ocean_v, ocean_f);
	generate_sphere(sphere_v, sphere_f);

	GLuint cube_vao = make_big_cube();
	GLuint cube_map_texture;
	create_cube_map(FRONT, BACK, TOP, BOTTOM, LEFT, RIGHT, &cube_map_texture );

	glm::vec4 min_bounds = glm::vec4(std::numeric_limits<float>::max());
	glm::vec4 max_bounds = glm::vec4(-std::numeric_limits<float>::max());
	for (const auto& vert : obj_vertices) {
		min_bounds = glm::min(vert, min_bounds);
		max_bounds = glm::max(vert, max_bounds);
	}
	std::cout << "min_bounds = " << glm::to_string(min_bounds) << "\n";
	std::cout << "max_bounds = " << glm::to_string(max_bounds) << "\n";

	// Setup our VAO array.
	CHECK_GL_ERROR(glGenVertexArrays(kNumVaos, &g_array_objects[0]));

	// Switch to the VAO for Geometry.
	CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kGeometryVao]));

	// Generate buffer objects
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &g_buffer_objects[kGeometryVao][0]));

	// Setup vertex data in a VBO.
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kGeometryVao][kVertexBuffer]));
	// NOTE: We do not send anything right now, we just describe it to OpenGL.
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER,
				sizeof(float) * obj_vertices.size() * 4, obj_vertices.data(),
				GL_STATIC_DRAW));
	CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));

	// Setup element array buffer.
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kGeometryVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER,
				sizeof(uint32_t) * obj_faces.size() * 3,
				obj_faces.data(), GL_STATIC_DRAW));

	/*
 	 * By far, the geometry is loaded into g_buffer_objects[kGeometryVao][*].
	 * These buffers are binded to g_array_objects[kGeometryVao]
	 */

	// FIXME: load the floor into g_buffer_objects[kFloorVao][*],
	//        and bind these VBO to g_array_objects[kFloorVao]


    //Switch to the VAO for floor
	CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kFloorVao]));
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &g_buffer_objects[kFloorVao][0]));
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kFloorVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * floor_v.size() * 4, &floor_v[0], GL_STATIC_DRAW));
    CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kFloorVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * floor_f.size() * 3, &floor_f[0], GL_STATIC_DRAW));

    // Switch to the VAO for ocean.
	CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kOceanVao]));
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &g_buffer_objects[kOceanVao][0]));
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kOceanVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * ocean_v.size() * 4, &ocean_v[0], GL_STATIC_DRAW));
    CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kOceanVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * ocean_f.size() * 4, &ocean_f[0], GL_STATIC_DRAW));

	// Switch to the VAO for sphere.
	CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kSphereVao]));
	CHECK_GL_ERROR(glGenBuffers(kNumVbos, &g_buffer_objects[kSphereVao][0]));
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kSphereVao][kVertexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * sphere_v.size() * 4, &sphere_v[0], GL_STATIC_DRAW));
    CHECK_GL_ERROR(glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0));
	CHECK_GL_ERROR(glEnableVertexAttribArray(0));
	CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kSphereVao][kIndexBuffer]));
	CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * sphere_f.size() * 3, &sphere_f[0], GL_STATIC_DRAW));


	// Setup vertex shader.
	GLuint vertex_shader_id = 0;
	const char* vertex_source_pointer = vertex_shader;
	CHECK_GL_ERROR(vertex_shader_id = glCreateShader(GL_VERTEX_SHADER));
	CHECK_GL_ERROR(glShaderSource(vertex_shader_id, 1, &vertex_source_pointer, nullptr));
	glCompileShader(vertex_shader_id);
	CHECK_GL_SHADER_ERROR(vertex_shader_id);

	// Setup geometry shader.
	GLuint geometry_shader_id = 0;
	const char* geometry_source_pointer = geometry_shader;
	CHECK_GL_ERROR(geometry_shader_id = glCreateShader(GL_GEOMETRY_SHADER));
	CHECK_GL_ERROR(glShaderSource(geometry_shader_id, 1, &geometry_source_pointer, nullptr));
	glCompileShader(geometry_shader_id);
	CHECK_GL_SHADER_ERROR(geometry_shader_id);

	// Setup fragment shader.
	GLuint fragment_shader_id = 0;
	const char* fragment_source_pointer = fragment_shader;
	CHECK_GL_ERROR(fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(glShaderSource(fragment_shader_id, 1, &fragment_source_pointer, nullptr));
	glCompileShader(fragment_shader_id);
	CHECK_GL_SHADER_ERROR(fragment_shader_id);

	// Let's create our program.
	GLuint program_id = 0;
	CHECK_GL_ERROR(program_id = glCreateProgram());
	CHECK_GL_ERROR(glAttachShader(program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(program_id, fragment_shader_id));
	CHECK_GL_ERROR(glAttachShader(program_id, geometry_shader_id));

	// Bind attributes.
	CHECK_GL_ERROR(glBindAttribLocation(program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(program_id, 0, "fragment_color"));
	glLinkProgram(program_id);
	CHECK_GL_PROGRAM_ERROR(program_id);

	// Get the uniform locations.
	GLint projection_matrix_location = 0;
	CHECK_GL_ERROR(projection_matrix_location =
			glGetUniformLocation(program_id, "projection"));
	GLint view_matrix_location = 0;
	CHECK_GL_ERROR(view_matrix_location =
			glGetUniformLocation(program_id, "view"));
	GLint light_position_location = 0;
	CHECK_GL_ERROR(light_position_location =
			glGetUniformLocation(program_id, "light_position"));

	// Setup fragment shader for the floor
	GLuint floor_fragment_shader_id = 0;
	const char* floor_fragment_source_pointer = floor_fragment_shader;
	CHECK_GL_ERROR(floor_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(glShaderSource(floor_fragment_shader_id, 1,
				&floor_fragment_source_pointer, nullptr));
	glCompileShader(floor_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(floor_fragment_shader_id);

	//Setup tessellation control shader for the floor
	GLuint floor_tess_control_shader_id = 0;
	const char* floor_tess_control_source_pointer = floor_tess_control_shader;
	CHECK_GL_ERROR(floor_tess_control_shader_id = glCreateShader(GL_TESS_CONTROL_SHADER));
	CHECK_GL_ERROR(glShaderSource(floor_tess_control_shader_id, 1,
				&floor_tess_control_source_pointer, nullptr));
	glCompileShader(floor_tess_control_shader_id);
	CHECK_GL_SHADER_ERROR(floor_tess_control_shader_id);

	//Setup tessellation evaluation shader for the floor
	GLuint floor_tess_evaluation_shader_id = 0;
	const char* floor_tess_evaluation_source_pointer = floor_tess_evaluation_shader;
	CHECK_GL_ERROR(floor_tess_evaluation_shader_id = glCreateShader(GL_TESS_EVALUATION_SHADER));
	CHECK_GL_ERROR(glShaderSource(floor_tess_evaluation_shader_id, 1,
				&floor_tess_evaluation_source_pointer, nullptr));
	glCompileShader(floor_tess_evaluation_shader_id);
	CHECK_GL_SHADER_ERROR(floor_tess_evaluation_shader_id);



	// Setup fragment shader for the ocean
	GLuint ocean_fragment_shader_id = 0;
	const char* ocean_fragment_source_pointer = ocean_fragment_shader;
	CHECK_GL_ERROR(ocean_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(glShaderSource(ocean_fragment_shader_id, 1,
				&ocean_fragment_source_pointer, nullptr));
	glCompileShader(ocean_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(ocean_fragment_shader_id);

	//Setup tessellation control shader for the ocean
	GLuint ocean_tess_control_shader_id = 0;
	const char* ocean_tess_control_source_pointer = ocean_tess_control_shader;
	CHECK_GL_ERROR(ocean_tess_control_shader_id = glCreateShader(GL_TESS_CONTROL_SHADER));
	CHECK_GL_ERROR(glShaderSource(ocean_tess_control_shader_id, 1,
				&ocean_tess_control_source_pointer, nullptr));
	glCompileShader(ocean_tess_control_shader_id);
	CHECK_GL_SHADER_ERROR(ocean_tess_control_shader_id);

	//Setup tessellation evaluation shader for the ocean
	GLuint ocean_tess_evaluation_shader_id = 0;
	const char* ocean_tess_evaluation_source_pointer = ocean_tess_evaluation_shader;
	CHECK_GL_ERROR(ocean_tess_evaluation_shader_id = glCreateShader(GL_TESS_EVALUATION_SHADER));
	CHECK_GL_ERROR(glShaderSource(ocean_tess_evaluation_shader_id, 1,
				&ocean_tess_evaluation_source_pointer, nullptr));
	glCompileShader(ocean_tess_evaluation_shader_id);
	CHECK_GL_SHADER_ERROR(ocean_tess_evaluation_shader_id);

	//Setup geometry shader for the ocean
	GLuint ocean_geometry_shader_id = 0;
	const char* ocean_geometry_source_pointer = ocean_geometry_shader;
	CHECK_GL_ERROR(ocean_geometry_shader_id = glCreateShader(GL_GEOMETRY_SHADER));
	CHECK_GL_ERROR(glShaderSource(ocean_geometry_shader_id, 1,
				&ocean_geometry_source_pointer, nullptr));
	glCompileShader(ocean_geometry_shader_id);
	CHECK_GL_SHADER_ERROR(ocean_geometry_shader_id);


	//Setup sphere fragment shader
	GLuint sphere_fragment_shader_id = 0;
	const char* sphere_fragment_source_pointer = sphere_fragment_shader;
	CHECK_GL_ERROR(sphere_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(glShaderSource(sphere_fragment_shader_id, 1, &sphere_fragment_source_pointer, nullptr));
	glCompileShader(sphere_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(sphere_fragment_shader_id);



	// Setup skybox vertex shader.
	GLuint skybox_vertex_shader_id = 0;
	const char* skybox_vertex_source_pointer = skybox_vertex_shader;
	CHECK_GL_ERROR(skybox_vertex_shader_id = glCreateShader(GL_VERTEX_SHADER));
	CHECK_GL_ERROR(glShaderSource(skybox_vertex_shader_id, 1, &skybox_vertex_source_pointer, nullptr));
	glCompileShader(skybox_vertex_shader_id);
	CHECK_GL_SHADER_ERROR(skybox_vertex_shader_id);

	// Setup skybox fragment shader.
	GLuint skybox_fragment_shader_id = 0;
	const char* skybox_fragment_source_pointer = skybox_fragment_shader;
	CHECK_GL_ERROR(skybox_fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER));
	CHECK_GL_ERROR(glShaderSource(skybox_fragment_shader_id, 1, &skybox_fragment_source_pointer, nullptr));
	glCompileShader(skybox_fragment_shader_id);
	CHECK_GL_SHADER_ERROR(skybox_fragment_shader_id);



	// FIXME: Setup another program for the floor, and get its locations.
	// Note: you can reuse the vertex and geometry shader objects
	GLuint floor_program_id = glCreateProgram();
	GLint floor_projection_matrix_location = 0;
	GLint floor_view_matrix_location = 0;
	GLint floor_light_position_location = 0;
	GLint floor_tess_outer_location = 0;
	GLint floor_tess_inner_location = 0;
	GLint floor_if_frame_location = 0;

    CHECK_GL_ERROR(floor_program_id = glCreateProgram());
	CHECK_GL_ERROR(glAttachShader(floor_program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(floor_program_id, floor_fragment_shader_id));
	CHECK_GL_ERROR(glAttachShader(floor_program_id, floor_tess_control_shader_id));
	CHECK_GL_ERROR(glAttachShader(floor_program_id, floor_tess_evaluation_shader_id));
	CHECK_GL_ERROR(glAttachShader(floor_program_id, geometry_shader_id));
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kFloorVao][kVertexBuffer]));
	//CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * floor_v.size() * 4, nullptr, GL_STATIC_DRAW));
    CHECK_GL_ERROR(glBindAttribLocation(floor_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(floor_program_id, 0, "fragment_color"));
	glLinkProgram(floor_program_id);
	CHECK_GL_PROGRAM_ERROR(floor_program_id);
	CHECK_GL_ERROR(floor_projection_matrix_location = glGetUniformLocation(floor_program_id, "projection"));
	CHECK_GL_ERROR(floor_view_matrix_location = glGetUniformLocation(floor_program_id, "view"));
	CHECK_GL_ERROR(floor_light_position_location = glGetUniformLocation(floor_program_id, "light_position"));
	CHECK_GL_ERROR(floor_tess_outer_location = glGetUniformLocation(floor_program_id, "tessellation_level_outer"));
	CHECK_GL_ERROR(floor_tess_inner_location = glGetUniformLocation(floor_program_id, "tessellation_level_inner"));
    CHECK_GL_ERROR(floor_if_frame_location = glGetUniformLocation(floor_program_id, "if_render_frame"));

    // for ocean program
	GLuint ocean_program_id = glCreateProgram();
	GLint ocean_projection_matrix_location = 0;
	GLint ocean_view_matrix_location = 0;
	GLint ocean_light_position_location = 0;
	GLint ocean_tess_outer_location = 0;
	GLint ocean_tess_inner_location = 0;
	GLint ocean_if_frame_location = 0;
	GLint ocean_time_location = 0;
	GLint ocean_tidal_time_location = 0;
	GLint ocean_if_tidal_location = 0;

    CHECK_GL_ERROR(ocean_program_id = glCreateProgram());
	CHECK_GL_ERROR(glAttachShader(ocean_program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(ocean_program_id, ocean_fragment_shader_id));
	CHECK_GL_ERROR(glAttachShader(ocean_program_id, ocean_tess_control_shader_id));
	CHECK_GL_ERROR(glAttachShader(ocean_program_id, ocean_tess_evaluation_shader_id));
	CHECK_GL_ERROR(glAttachShader(ocean_program_id, ocean_geometry_shader_id));
	CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kOceanVao][kVertexBuffer]));
	//CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * floor_v.size() * 4, nullptr, GL_STATIC_DRAW));
    CHECK_GL_ERROR(glBindAttribLocation(ocean_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(ocean_program_id, 0, "fragment_color"));
	glLinkProgram(ocean_program_id);
	CHECK_GL_PROGRAM_ERROR(ocean_program_id);
	CHECK_GL_ERROR(ocean_projection_matrix_location = glGetUniformLocation(ocean_program_id, "projection"));
	CHECK_GL_ERROR(ocean_view_matrix_location = glGetUniformLocation(ocean_program_id, "view"));
	CHECK_GL_ERROR(ocean_light_position_location = glGetUniformLocation(ocean_program_id, "light_position"));
	CHECK_GL_ERROR(ocean_tess_outer_location = glGetUniformLocation(ocean_program_id, "tessellation_level_outer"));
	CHECK_GL_ERROR(ocean_tess_inner_location = glGetUniformLocation(ocean_program_id, "tessellation_level_inner"));
	CHECK_GL_ERROR(ocean_if_frame_location = glGetUniformLocation(ocean_program_id, "if_render_frame"));
	CHECK_GL_ERROR(ocean_time_location = glGetUniformLocation(ocean_program_id, "time"));
	CHECK_GL_ERROR(ocean_tidal_time_location = glGetUniformLocation(ocean_program_id, "tidal_time"));
	CHECK_GL_ERROR(ocean_if_tidal_location = glGetUniformLocation(ocean_program_id, "if_generate_tidal"));


	//for sphere program


	GLuint sphere_program_id = 0;
	CHECK_GL_ERROR(sphere_program_id = glCreateProgram());
	CHECK_GL_ERROR(glAttachShader(sphere_program_id, vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(sphere_program_id, sphere_fragment_shader_id));
	CHECK_GL_ERROR(glAttachShader(sphere_program_id, geometry_shader_id));



	// Bind attributes.
	CHECK_GL_ERROR(glBindAttribLocation(sphere_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(sphere_program_id, 0, "fragment_color"));
	glLinkProgram(sphere_program_id);
	CHECK_GL_PROGRAM_ERROR(sphere_program_id);

	// Get the uniform locations.
	GLint sphere_projection_matrix_location = 0;
	CHECK_GL_ERROR(sphere_projection_matrix_location =
			glGetUniformLocation(sphere_program_id, "projection"));
	GLint sphere_view_matrix_location = 0;
	CHECK_GL_ERROR(sphere_view_matrix_location =
			glGetUniformLocation(sphere_program_id, "view"));
	GLint sphere_light_position_location = 0;
	CHECK_GL_ERROR(sphere_light_position_location =
			glGetUniformLocation(sphere_program_id, "light_position"));


    // for skybox program
	GLuint skybox_program_id = 0;
	CHECK_GL_ERROR(skybox_program_id = glCreateProgram());
	CHECK_GL_ERROR(glAttachShader(skybox_program_id, skybox_vertex_shader_id));
	CHECK_GL_ERROR(glAttachShader(skybox_program_id, skybox_fragment_shader_id));


	// Bind attributes.
	CHECK_GL_ERROR(glBindAttribLocation(skybox_program_id, 0, "vertex_position"));
	CHECK_GL_ERROR(glBindFragDataLocation(skybox_program_id, 0, "fragment_color"));
	glLinkProgram(skybox_program_id);
	CHECK_GL_PROGRAM_ERROR(skybox_program_id);

	// Get the uniform locations.
	GLint skybox_projection_matrix_location = 0;
	CHECK_GL_ERROR(skybox_projection_matrix_location =
			glGetUniformLocation(skybox_program_id, "projection"));
	GLint skybox_view_matrix_location = 0;
	CHECK_GL_ERROR(skybox_view_matrix_location =
			glGetUniformLocation(skybox_program_id, "view"));



	struct timespec start_time;
	clock_gettime(CLOCK_REALTIME, &start_time);

	glm::vec4 light_position = glm::vec4(-10.0f, 10.0f, 0.0f, 1.0f);
	float aspect = 0.0f;
	float theta = 0.0f;

	struct timespec tidal_start_time;
	clock_gettime(CLOCK_REALTIME, &tidal_start_time);

	while (!glfwWindowShouldClose(window)) {
		// Setup some basic window stuff.
		glfwGetFramebufferSize(window, &window_width, &window_height);
		glViewport(0, 0, window_width, window_height);
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glEnable(GL_DEPTH_TEST);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glDepthFunc(GL_LESS);

		if (if_save) {
			SaveObj("geometry.obj", obj_vertices, obj_faces);
			if_save = false;
		}



		// Switch to the Geometry VAO.
		CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kGeometryVao]));

		if (g_menger && g_menger->is_dirty()) {
			obj_vertices.clear();
			obj_faces.clear();
			g_menger->generate_geometry(obj_vertices, obj_faces);
			g_menger->set_clean();

			// FIXME: Upload your vertex data here.
			CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kGeometryVao][kVertexBuffer]));
			CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * obj_vertices.size() * 4, &obj_vertices[0], GL_STATIC_DRAW));
		}

		// Compute the projection matrix.
		aspect = static_cast<float>(window_width) / window_height;
		glm::mat4 projection_matrix =
			glm::perspective(glm::radians(45.0f), aspect, 0.0001f, 1000.0f);

		// Compute the view matrix
		// FIXME: change eye and center through mouse/keyboard events.
		glm::mat4 view_matrix = g_camera.get_view_matrix();

		// Use our program.
		CHECK_GL_ERROR(glUseProgram(program_id));

		// Pass uniforms in.
		CHECK_GL_ERROR(glUniformMatrix4fv(projection_matrix_location, 1, GL_FALSE,
					&projection_matrix[0][0]));
		CHECK_GL_ERROR(glUniformMatrix4fv(view_matrix_location, 1, GL_FALSE,
					&view_matrix[0][0]));
		CHECK_GL_ERROR(glUniform4fv(light_position_location, 1, &light_position[0]));

		// Draw our triangles.
		CHECK_GL_ERROR(glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)); 
		CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, obj_faces.size() * 3, GL_UNSIGNED_INT, 0));

		// FIXME: Render the floor
		// Note: What you need to do is
		// 	1. Switch VAO
		// 	2. Switch Program
		// 	3. Pass Uniforms
		// 	4. Call glDrawElements, since input geometry is
		// 	indicated by VAO.


		if (if_faces_render){
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		} else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}

	    //switch to Floor VAO
		CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kFloorVao]));
        //bind floor vetices
		CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kFloorVao][kVertexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * floor_v.size() * 4, &floor_v[0], GL_STATIC_DRAW));
		//bind floor faces
		CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kFloorVao][kIndexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * floor_f.size()*3, &floor_f[0], GL_STATIC_DRAW));
		//switch program
		CHECK_GL_ERROR(glUseProgram(floor_program_id));

		//pass uniforms
		CHECK_GL_ERROR(glUniformMatrix4fv(floor_projection_matrix_location, 1, GL_FALSE, &projection_matrix[0][0]));
		CHECK_GL_ERROR(glUniformMatrix4fv(floor_view_matrix_location, 1, GL_FALSE, &view_matrix[0][0]));
		CHECK_GL_ERROR(glUniform4fv(floor_light_position_location, 1, &light_position[0]));
		CHECK_GL_ERROR(glUniform1f(floor_tess_outer_location, tessellation_level_outer));
		CHECK_GL_ERROR(glUniform1f(floor_tess_inner_location, tessellation_level_inner));
		CHECK_GL_ERROR(glUniform1i(floor_if_frame_location, if_render_frame));

        //draw elements of floor
        CHECK_GL_ERROR(glPatchParameteri(GL_PATCH_VERTICES, 3));
        if (!if_ocean_mode){
        	CHECK_GL_ERROR(glDrawElements(GL_PATCHES, floor_f.size() * 3, GL_UNSIGNED_INT, 0));
        }


	    //switch to Ocean VAO
		CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kOceanVao]));
        //bind floor vetices
		CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kOceanVao][kVertexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * ocean_v.size() * 4, &ocean_v[0], GL_STATIC_DRAW));
		//bind floor faces
		CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kOceanVao][kIndexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * ocean_f.size()*4, &ocean_f[0], GL_STATIC_DRAW));
		//switch program
		CHECK_GL_ERROR(glUseProgram(ocean_program_id));

		//pass uniforms
		CHECK_GL_ERROR(glUniformMatrix4fv(ocean_projection_matrix_location, 1, GL_FALSE, &projection_matrix[0][0]));
		CHECK_GL_ERROR(glUniformMatrix4fv(ocean_view_matrix_location, 1, GL_FALSE, &view_matrix[0][0]));
		CHECK_GL_ERROR(glUniform4fv(ocean_light_position_location, 1, &light_position[0]));
		CHECK_GL_ERROR(glUniform1f(ocean_tess_outer_location, tessellation_level_outer));
		CHECK_GL_ERROR(glUniform1f(ocean_tess_inner_location, tessellation_level_inner));
		CHECK_GL_ERROR(glUniform1i(ocean_if_frame_location, if_render_frame));
		CHECK_GL_ERROR(glUniform1i(ocean_if_tidal_location, if_generate_tidal));

        if (if_generate_tidal && !if_reset_tidal_time){
        	//if_generate_tidal = false;
        	if_reset_tidal_time = true;
        	clock_gettime(CLOCK_REALTIME, &tidal_start_time);
        }

        struct timespec current_time;
        clock_gettime(CLOCK_REALTIME, &current_time);
        float time_diff = (current_time.tv_sec - start_time.tv_sec) + (float(current_time.tv_nsec - start_time.tv_nsec))/(1000000000LL);
        float tidal_time_diff = (current_time.tv_sec - tidal_start_time.tv_sec) + (float(current_time.tv_nsec - tidal_start_time.tv_nsec))/(1000000000LL);
        CHECK_GL_ERROR(glUniform1f(ocean_time_location, time_diff));
        CHECK_GL_ERROR(glUniform1f(ocean_tidal_time_location, tidal_time_diff));
        //draw elements of ocean
        CHECK_GL_ERROR(glPatchParameteri(GL_PATCH_VERTICES, 4));
        if (if_ocean_mode) {
			CHECK_GL_ERROR(glDrawElements(GL_PATCHES, ocean_f.size() * 4, GL_UNSIGNED_INT, 0));
        }

        if (tidal_time_diff > 60){
        	if_generate_tidal = false;
        }


        //switch to Sphere VAO
		CHECK_GL_ERROR(glBindVertexArray(g_array_objects[kSphereVao]));
        //bind floor vetices
		CHECK_GL_ERROR(glBindBuffer(GL_ARRAY_BUFFER, g_buffer_objects[kSphereVao][kVertexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * sphere_v.size() * 4, &sphere_v[0], GL_STATIC_DRAW));
		//bind floor faces
		CHECK_GL_ERROR(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_buffer_objects[kSphereVao][kIndexBuffer]));
		CHECK_GL_ERROR(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * sphere_f.size()*3, &sphere_f[0], GL_STATIC_DRAW));
		//switch program
		CHECK_GL_ERROR(glUseProgram(sphere_program_id));

		//pass uniforms
		CHECK_GL_ERROR(glUniformMatrix4fv(sphere_projection_matrix_location, 1, GL_FALSE, &projection_matrix[0][0]));
		CHECK_GL_ERROR(glUniformMatrix4fv(sphere_view_matrix_location, 1, GL_FALSE, &view_matrix[0][0]));
		CHECK_GL_ERROR(glUniform4fv(sphere_light_position_location, 1, &light_position[0]));

        //draw elements of sphere
        //CHECK_GL_ERROR(glPatchParameteri(GL_PATCH_VERTICES, 4));
        CHECK_GL_ERROR(glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)); 
        CHECK_GL_ERROR(glDrawElements(GL_TRIANGLES, sphere_f.size() * 3, GL_UNSIGNED_INT, 0));


        // for cubemap program
        CHECK_GL_ERROR(glUseProgram(skybox_program_id));
        CHECK_GL_ERROR(glUniformMatrix4fv(skybox_projection_matrix_location, 1, GL_FALSE, &projection_matrix[0][0]));
		CHECK_GL_ERROR(glUniformMatrix4fv(skybox_view_matrix_location, 1, GL_FALSE, &view_matrix[0][0]));
		glActiveTexture( GL_TEXTURE0 );
		glBindTexture( GL_TEXTURE_CUBE_MAP, cube_map_texture );
		glBindVertexArray( cube_vao );
		glDrawArrays( GL_TRIANGLES, 0, 36 );
        glDepthMask( GL_TRUE );



		// Poll and swap.
		glfwPollEvents();
		glfwSwapBuffers(window);
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}




