
#ifndef WEBGL_NOISE_COMMON_GLSL
#define WEBGL_NOISE_COMMON_GLSL

vec4 mod289(vec4 x)
{
	return x - floor(x * (float(1.0) / float(289.0))) * float(289.0);
}

// Modulo 289 without a division (only multiplications)
vec3 mod289(vec3 x)
{
	return x - floor(x * (float(1.0) / float(289.0))) * float(289.0);
}

vec2 mod289(vec2 x)
{
	return x - floor(x * (float(1.0) / float(289.0))) * float(289.0);
}

float mod289(float x)
{
	return x - floor(x * (float(1.0) / float(289.0))) * float(289.0);
}

vec4 permute(vec4 x) { return mod289(((x * float(34.0)) + float(10.0)) * x); }

// Permutation polynomial: (34x^2 + 10x) mod 289
vec3 permute(vec3 x) { return mod289((float(34.0) * x + float(10.0)) * x); }

vec2 permute(vec2 x) { return mod289(((x * float(34.0)) + float(10.0)) * x); }

float permute(float x) { return mod289(((x * float(34.0)) + float(10.0)) * x); }

vec4 taylorInvSqrt(vec4 r)
{
	return float(1.79284291400159) - float(0.85373472095314) * r;
}

float taylorInvSqrt(float r)
{
	return float(1.79284291400159) - float(0.85373472095314) * r;
}

// Modulo 7 without a division
vec4 mod7(vec4 x)
{
	return x - floor(x * (float(1.0) / float(7.0))) * float(7.0);
}

// Modulo 7 without a division
vec3 mod7(vec3 x)
{
	return x - floor(x * (float(1.0) / float(7.0))) * float(7.0);
}

// Modulo 7 without a division
vec2 mod7(vec2 x)
{
	return x - floor(x * (float(1.0) / float(7.0))) * float(7.0);
}

// Modulo 7 without a division
float mod7(float x)
{
	return x - floor(x * (float(1.0) / float(7.0))) * float(7.0);
}

#endif
