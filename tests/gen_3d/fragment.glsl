#version 430 core

uniform vec3 scale;

in vec4 out_color;
in vec3 pos;
in vec3 intPos;
in vec3 normal;
in vec3 trueNormal;
in vec3 triangleVert0Pos;

out vec4 FragColor;

const float factorGrid = 0.0;

void main()
{
	vec3 f = fract(intPos);
	if ((f.x < factorGrid) || (f.z < factorGrid)) {
		FragColor = vec4(1, 0, 0, 1);
	} else {

		if (scale.y * 0.07 > triangleVert0Pos.y) {
			FragColor = vec4(0.1, 0.1, out_color.x * 0.3 + 0.5, 1);
		} else if (normal.y < 0.98) {
			FragColor = vec4(out_color.xyz * 0.3 + 0.4, 1);
		} else {
			if (scale.y * 0.085 > triangleVert0Pos.y) {
				FragColor = vec4(out_color.xyz * 0.2 + 0.8, 1);
				FragColor.b = 0.1;
			} else if (scale.y * 0.62 < triangleVert0Pos.y) {
				FragColor = vec4(out_color.xyz * 0.2 + 0.8, 1);
			} else {
				FragColor = vec4(0.1, out_color.x * 0.3 + 0.7, 0.1, 1);
			}
		}

		float light = dot(trueNormal, normalize(vec3(-0.6, -0.1, -0.5)));
// 		light *= light * light;
		light = light * 0.5 + 0.5;
		light = light * light;
		light = light * 1.2 - 0.2;
		light = 1.0 - light*0.6;

		FragColor = vec4(FragColor.xyz * light, 1);
	}
}
