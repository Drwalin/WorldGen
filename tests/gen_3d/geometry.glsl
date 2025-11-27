#version 430 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

uniform vec3 scale;
uniform ivec2 size;

in vec4 _in_out_color[3];
in vec3 _in_pos[3];
in float _in_water[3];
in vec2 _in_uv[3];

out vec4 out_color;
out vec3 lightNormal;
out vec3 pos;
out vec2 uv;
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
		uv = _in_uv[i];
		out_color = _in_out_color[i];
		outWater = _in_water[i];
		lightNormal = n;
		EmitVertex();
	}
	EndPrimitive();
}
