#ifndef C_GLSL_RANDOM_H
#define C_GLSL_RANDOM_H

#include "../include/worldgen/GlslSharedWithCpp_Header.h"

#if GL_core_profile
#else
using uint = unsigned int;
namespace glsl_noise
{
using namespace glm;
#endif

// Random function code based on:
//      https://www.reedbeta.com/blog/hash-functions-for-gpu-rendering/

inline uint pcg_hash(uint seed)
{
	uint state = seed * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}

inline uint pcg_hash_seed(uint rng_state, uint seed)
{
	return (rng_state + seed) * 747796405u;
}

inline uvec2 rand_pcg(uint rng_state)
{
	uint state = rng_state;
	rng_state = rng_state * 747796405u + 2891336453u;
	uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return uvec2((word >> 22u) ^ word, rng_state);
}

inline uvec4 RandomUint(ivec3 p)
{
	uint state = pcg_hash_seed(pcg_hash_seed(pcg_hash(p.x), p.y), p.z);
	uvec4 ret;
	uvec2 h;
	h = rand_pcg(state);
	ret.x = h.x;
	h = rand_pcg(h.y);
	ret.y = h.x;
	h = rand_pcg(h.y);
	ret.z = h.x;
	h = rand_pcg(h.y);
	ret.w = h.x;
	return ret;
}

inline ivec4 RandomInt(ivec3 p) { return ivec4(RandomUint(p)); }

#if GL_core_profile
#else
}
#endif

#endif
