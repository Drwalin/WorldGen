#version 430 core

uniform sampler2D colorTex;
uniform sampler2D heightTex;

uniform ivec2 size;
uniform vec3 scale;
uniform int useWater = 0;
uniform int gridWidth;

// in vec4 out_color;
in vec3 pos;
in vec2 uv;
in vec3 lightNormal;
in vec3 triangleVert0Pos;
// in float outWater;

out vec4 FragColor;

const float _heightScale = 350.0;
uniform float generatorYScale = 1.0;


vec3 CalcNormal(vec2 coords) {
	float dx = 2.0 / float(size.x);
	vec3 p[3];
	p[0].xz = coords + vec2(dx,dx);
	p[1].xz = coords + vec2(dx,-dx);
	p[2].xz = coords + vec2(-dx,0);
	
	for (int i=0; i<3; ++i) {
		p[i].y = texture(heightTex, p[i].xz).x;
		p[i].xz *= vec2(size);
		p[i] = p[i] * scale;
	}
	
	vec3 n = normalize(cross(p[1] - p[0], p[2] - p[0]));
	if (n.y < 0.0) {
		n = -n;
	}
	return n;
}

void main()
{
	float heightScale = _heightScale * generatorYScale;
	
	vec2 hw = texture(heightTex, uv).xy;
	vec3 normal = CalcNormal(uv);
	vec4 out_color = texture(colorTex, uv);
	
	float height = hw.x * scale.y;
	float water = hw.y;
	vec4 colorSand = vec4(out_color.xyz * 0.2 + 0.8, 1);
	colorSand.b = 0.1;
	vec4 colorStone = vec4(out_color.xyz * 0.3 + 0.4, 1);
	vec4 colorSnow = vec4(out_color.xyz * 0.2 + 0.8, 1);
	vec4 colorGrass = vec4(0.1, out_color.x * 0.3 + 0.7, 0.1, 1);
	vec4 colorWater = vec4(0.1, 0.1, out_color.x * 0.3 + 0.5, 1);
	
	/*
	if (scale.y * 0.07 > height) {
		// water
		FragColor = vec4(0.1, 0.1, out_color.x * 0.3 + 0.5, 1);
	} else
	*/
	if (normal.y < 0.9) {
		// stone
		FragColor = colorStone;
	} else {
		if (scale.y * 0.085 * heightScale > height) {
			// sand
			FragColor = colorSand;
		} else if (scale.y * 0.62 * heightScale < height) {
			// snow
			FragColor = colorSnow;
		} else {
			// grass
			FragColor = colorGrass;
		}
	}
	
	
	if (useWater == 1) {
		float w = clamp(water, 0.0, 1.0);
		if (w < 0.25) {
			w *= 10.0;
			w = sqrt(w);
			w = sqrt(w);
			w *= 0.1;
		}
		if (w < 0.05) {
			discard;
		}
		FragColor = mix(FragColor, colorWater, w);
	} else {
		float w = clamp(water, 0.0, 1.0);
		w = sqrt(w);
// 		if (w < 0.25) {
// 			w *= 10.0;
// 			w = sqrt(w);
// 			w *= 0.1;
// 		}
		FragColor = mix(FragColor, colorSand, w);
	}

	float light = dot(lightNormal, normalize(vec3(-0.6, -0.1, -0.5)));
// 	light *= light * light;
	
	light = light * 0.5 + 0.5;
// 	light = light * light;
// 	light = light * 1.2 - 0.2;
// 	light = 1.0 - light*0.6;

	FragColor = vec4(FragColor.xyz * light, 1.0 - float(useWater) * 0.3);
	if (gridWidth > 1600) {
		FragColor.w = 1.0;
	}
}
