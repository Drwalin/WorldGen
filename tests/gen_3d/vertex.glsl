#version 430 core

uniform sampler2D colorTex;
uniform sampler2D heightTex;

uniform ivec2 meshSize;

uniform vec3 cameraPos;
uniform vec2 centerOffset;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 renderScaleVec = vec3(1,1,1);
uniform int useWater = 0;

uniform ivec2 size;
uniform vec3 scale;

out vec4 _in_out_color;
out vec3 _in_pos;
out float _in_water;
out vec2 _in_uv;


float f0(float f) {
	return -f*f*0.5 + 1.5 * f;
}

float f1(float f, float d) {
	return f*f*(d-1.0) + f * (2.0 - d);
}

void main()
{
	int _x = gl_VertexID % meshSize.x;
	int _z = gl_VertexID / meshSize.x;
	
	vec2 coord = vec2(_x, _z) / (meshSize-1);
	
	vec2 meshScale = vec2(meshSize) / vec2(size);
	
	vec2 pos;
	vec2 minBorder = centerOffset - meshScale * 0.25;
	vec2 maxBorder = centerOffset + meshScale * 0.25;
	
	pos = (coord - 0.25) * meshScale + minBorder;
	
	if (coord.x < 0.25) {
		float d = meshScale.x * 0.25 / minBorder.x;
		float f = coord.x * 4.0;
		f = f1(f, d);
		pos.x = f * minBorder.x;
	} else if (coord.x > 0.75) {
		float d = meshScale.x * 0.25 / (1.0 - maxBorder.x);
		float f = (coord.x - 0.75) * 4.0;
		f = 1.0 - f1(1.0 - f, d);
		pos.x = f * (1.0f - maxBorder.x) + maxBorder.x; 
	}
	
	if (coord.y < 0.25) {
		float d = meshScale.y * 0.25 / minBorder.y;
		float f = coord.y * 4.0;
		f = f1(f, d);
		pos.y = f * minBorder.y;
	} else if (coord.y > 0.75) {
		float d = meshScale.y * 0.25 / (1.0 - maxBorder.y);
		float f = (coord.y - 0.75) * 4.0;
		f = 1.0 - f1(1.0 - f, d);
		pos.y = f * (1.0f - maxBorder.y) + maxBorder.y; 
	}
	
	pos = clamp(pos, 0.0, 1.0);
	
	float x = int(pos.x * float(size.x));
	float z = int(pos.y * float(size.y));
	
	pos = vec2(x, z) / size;
	coord = pos;
	
	vec2 heights = texture(heightTex, coord).xy;
	float height = heights.x;
	vec4 color = texture(colorTex, coord);
	float water = heights.y;
	
	_in_uv = coord;
	_in_pos = vec3(x * scale.x, (height + water * float(useWater)) * scale.y, z * scale.z);
	gl_Position = projection * view * model * vec4(_in_pos * renderScaleVec, 1);
	_in_out_color = color;
	_in_water = water;
}

