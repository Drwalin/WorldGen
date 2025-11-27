#version 430 core

uniform vec3 scale;
uniform int useWater = 0;
uniform int gridWidth;

in vec4 out_color;
in vec3 pos;
in vec3 intPos;
in vec3 normal;
in vec3 triangleVert0Pos;
in float outWater;

out vec4 FragColor;

const float factorGrid = 0.0;
const float heightScale = 350.0;

void main()
{
	vec4 colorSand = vec4(out_color.xyz * 0.2 + 0.8, 1);
	colorSand.b = 0.1;
	vec4 colorStone = vec4(out_color.xyz * 0.3 + 0.4, 1);
	vec4 colorSnow = vec4(out_color.xyz * 0.2 + 0.8, 1);
	vec4 colorGrass = vec4(0.1, out_color.x * 0.3 + 0.7, 0.1, 1);
	vec4 colorWater = vec4(0.1, 0.1, out_color.x * 0.3 + 0.5, 1);
	
	
	
	vec3 f = fract(intPos);
	if ((f.x < factorGrid) || (f.z < factorGrid)) {
		FragColor = vec4(1, 0, 0, 1);
	} else {
		/*
		if (scale.y * 0.07 > triangleVert0Pos.y) {
			// water
			FragColor = vec4(0.1, 0.1, out_color.x * 0.3 + 0.5, 1);
		} else
		*/
		if (normal.y < 0.9) {
			// stone
			FragColor = colorStone;
		} else {
			if (scale.y * 0.085 * heightScale > triangleVert0Pos.y) {
				// sand
				FragColor = colorSand;
			} else if (scale.y * 0.62 * heightScale < triangleVert0Pos.y) {
				// snow
				FragColor = colorSnow;
			} else {
				// grass
				FragColor = colorGrass;
			}
		}
		
		
		if (useWater == 1) {
			float w = clamp(outWater, 0.0, 1.0);
			if (w < 0.25) {
				w *= 10.0;
				w = sqrt(w);
				w = sqrt(w);
				w *= 0.1;
			}
			if (gridWidth <= 1600) {
				if (w < 0.05) {
					discard;
				}
			}
			FragColor = mix(FragColor, colorWater, w);
		} else {
			float w = clamp(outWater, 0.0, 1.0);
			w = sqrt(w);
// 			if (w < 0.25) {
// 				w *= 10.0;
// 				w = sqrt(w);
// 				w *= 0.1;
// 			}
			FragColor = mix(FragColor, colorSand, w);
		}

		float light = dot(normal, normalize(vec3(-0.6, -0.1, -0.5)));
// 		light *= light * light;
		
		light = light * 0.5 + 0.5;
// 		light = light * light;
// 		light = light * 1.2 - 0.2;
// 		light = 1.0 - light*0.6;

		FragColor = vec4(FragColor.xyz * light, 1.0 - float(useWater) * 0.3);
		if (gridWidth > 1600) {
			FragColor.w = 1.0;
		}
	}
}
