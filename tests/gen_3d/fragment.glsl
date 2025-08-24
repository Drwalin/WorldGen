#version 430 core

in vec4 out_color;

out vec4 FragColor;

void main() {
	FragColor = vec4(1,1,1,1);
	FragColor *= vec4(out_color.xyz, 1);
}

