#version 430 core

in vec4 out_color;
in vec3 pos;
in vec3 intPos;

out vec4 FragColor;

const float factorGrid = 0.0;

void main() {
	vec3 f = fract(intPos);
	if ((f.x < factorGrid) || (f.z < factorGrid)) {
		FragColor = vec4(1,0,0,1);
	} else {
		FragColor = vec4(out_color.xyz, 1);
	}
}

