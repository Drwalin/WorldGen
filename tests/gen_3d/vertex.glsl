#version 430 core

layout ( location = 0 ) in float height;
layout ( location = 1 ) in vec4 color;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform ivec2 size;
uniform vec3 scale;

out vec2 out_uv;
out vec4 out_color;

void main() {
	int x = gl_VertexID % size.x;
	int z = gl_VertexID / size.x;
	gl_Position = projection * view * model * vec4(vec3(x, height, z)*scale, 1);
	out_color = color;
}

