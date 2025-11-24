#version 430 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

uniform vec3 scale;

in vec4 _in_out_color[3];
in vec3 _in_pos[3];
in vec3 _in_intPos[3];
in float _in_water[3];

out vec4 out_color;
out vec3 intPos;
out vec3 normal;
out vec3 pos;
out vec3 triangleVert0Pos;
out float outWater;

void main()
{
	vec3 a = _in_pos[0];
	vec3 b = _in_pos[1];
	vec3 c = _in_pos[2];
	triangleVert0Pos = _in_pos[0];

	vec3 n = normalize(cross(b - a, c - a));
	if (n.y < 0) {
		n = -n;
	}

	for (int i = 0; i < 3; ++i) {
		gl_Position = gl_in[i].gl_Position;
		pos = _in_pos[i];
		out_color = _in_out_color[i];
		intPos = _in_intPos[i];
		normal = n;
		outWater = _in_water[i];
		EmitVertex();
	}
	EndPrimitive();
}
