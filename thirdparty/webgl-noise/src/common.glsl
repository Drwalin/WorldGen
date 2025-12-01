#ifndef WEBGL_NOISE_COMMON_GLSL
#define WEBGL_NOISE_COMMON_GLSL

// Modulo 289 without a division (only multiplications)

#define DEFINE_MOD289(TYPE)                                                    \
	inline static TYPE mod289(TYPE x)                                          \
	{                                                                          \
		return x - floor(x * (float(1.0) / float(289.0))) * float(289.0);      \
	}
DEFINE_MOD289(vec4)
DEFINE_MOD289(vec3)
DEFINE_MOD289(vec2)
DEFINE_MOD289(float)
#undef DEFINE_MOD289

// Permutation polynomial: (34x^2 + 10x) mod 289
#define DEFINE_PERMUTE(TYPE)                                                   \
	inline static TYPE permute(TYPE x)                                         \
	{                                                                          \
		return mod289(((x * float(34.0)) + float(10.0)) * x);                  \
	}
DEFINE_PERMUTE(vec4)
DEFINE_PERMUTE(vec3)
DEFINE_PERMUTE(vec2)
DEFINE_PERMUTE(float)
#undef DEFINE_PERMUTE

#define DEFINE_TAYLORSQRT(TYPE)                                                \
	inline static TYPE taylorInvSqrt(TYPE r)                                   \
	{                                                                          \
		return float(1.79284291400159) - float(0.85373472095314) * r;          \
	}
DEFINE_TAYLORSQRT(vec4)
DEFINE_TAYLORSQRT(vec3)
DEFINE_TAYLORSQRT(vec2)
DEFINE_TAYLORSQRT(float)
#undef DEFINE_TAYLORSQRT

// Modulo 7 without a division
#define DEFINE_MOD7(TYPE)                                                      \
	inline static TYPE mod7(TYPE x)                                            \
	{                                                                          \
		return x - floor(x * (float(1.0) / float(7.0))) * float(7.0);          \
	}
DEFINE_MOD7(vec4)
DEFINE_MOD7(vec3)
DEFINE_MOD7(vec2)
DEFINE_MOD7(float)
#undef DEFINE_MOD7

#endif
