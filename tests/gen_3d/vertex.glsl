#version 430 core

layout(location = 0) in float height;
layout(location = 1) in vec4 color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform ivec2 size;
uniform vec3 scale;

out vec4 _in_out_color;
out vec3 _in_pos;
out vec3 _in_intPos;

void main()
{
	int x = gl_VertexID % size.x;
	int z = gl_VertexID / size.x;
	_in_intPos = vec3(x, height, z);
	_in_pos = vec3(x * scale.x, height, z * scale.z);
	gl_Position = projection * view * model * vec4(_in_pos, 1);
	_in_out_color = color;
}
