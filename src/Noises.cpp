// Based on Inigo Quilez articles and codes

// The MIT License
// Copyright (C) 2017 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright
// notice and this permission notice shall be included in all copies or
// substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
// TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// https://www.youtube.com/c/InigoQuilez
// https://iquilezles.org/

#include <random>

#include "../include/worldgen/Hash.hpp"
#include "../include/worldgen/Noises.hpp"

#include "../include/worldgen/GlslSharedWithCpp_Header.h"
#ifndef C_GLSL_SHARED_CODE_H
#endif

namespace glsl_noise
{
using namespace glm;
#include "../thirdparty/webgl-noise/src/classicnoise2D.glsl"
#include "../thirdparty/webgl-noise/src/cellular2D.glsl"
#include "../thirdparty/webgl-noise/src/noise4D.glsl"
#include "../thirdparty/psrdnoise/src/psrdnoise2.glsl"
#include "../thirdparty/psrdnoise/src/psrdnoise3.glsl"
} // namespace glsl_noise

namespace wg
{

/*
static glm::vec4 cubic(float v)
{
	glm::vec4 n = glm::vec4(1.0, 2.0, 3.0, 4.0) - v;
	glm::vec4 s = n * n * n;
	float x = s.x;
	float y = s.y - 4.0 * s.x;
	float z = s.z - 4.0 * s.y + 6.0 * s.x;
	float w = 6.0 - x - y - z;
	return glm::vec4(x, y, z, w) * (1.0f / 6.0f);
}
*/

static float sinc(float x)
{
	const float PI = 3.1415926;
	return sin(PI * x) / (PI * x);
}

static float LWeight(float x)
{
	if (fabs(x) < 1e-8) {
		return 1.0f;
	} else {
		return sinc(x) * sinc(x / 3);
	}
}

float Lanczos(glm::vec2 texCoord, float values[5][5])
{
	glm::vec2 imgCoord = texCoord;
	glm::vec2 samplePos = floor(imgCoord - 0.5f) + 0.5f;
	glm::vec2 interpFactor = imgCoord - samplePos;
	float nSum = 0;
	float fX, fY;
	for (int m = -2; m <= 3; m++) {
		fX = LWeight(float(m) - interpFactor.x);
		float vecCooef1 = fX;
		for (int n = -2; n <= 3; n++) {
			float vecData = values[m + 2][n + 2];
			fY = LWeight(float(n) - interpFactor.y);
			float vecCoeef2 = fY;
			nSum += vecData * vecCoeef2 * vecCooef1;
		}
	}
	return nSum;
}

float Lanczos(glm::vec2 texCoord, int seed)
{
	glm::vec2 imgCoord = texCoord;
	glm::vec2 samplePos = floor(imgCoord - 0.5f) + 0.5f;

	float values[5][5];
	for (int m = -2; m <= 3; m++) {
		for (int n = -2; n <= 3; n++) {
			values[m + 2][n + 2] = hx::Hash(samplePos + glm::vec2(m, n), seed);
		}
	}
	return Lanczos(texCoord, values);
}

static inline float CubicFactor(float f)
{
	float f2 = f * f;
	return -2 * f * f2 + 3 * f2;
}

/*
static inline float Cubic(float a, float b, float f)
{
	f = CubicFactor(f);
	return a * (1 - f) + b * f;
}
*/

SimplexNoise::SimplexNoise(uint64_t seed)
{
	Init(seed);
}

SimplexNoise::~SimplexNoise() {}

void SimplexNoise::Init(uint64_t seed)
{
	this->seed = seed;
	std::mt19937_64 mt(seed);
	auto gen = [&mt](float min, float max) -> float {
		return std::uniform_real_distribution<float>(min, max)(mt);
	};
	for (int i = 0; i < AMOUNT_SEEDED_VALS; ++i) {
		seedi[i] = mt();
		seedf[i] = gen(0, 1);
		seedrad[i] = gen(0, M_PI * 2);
		float rad = gen(0, M_PI * 2);
		float s = sin(rad);
		float c = cos(rad);
		seedrot[i] = glm::mat2(c, -s, s, c);
	}
}

float SimplexNoise::Noise(glm::vec2 p)
{
	glm::vec2 grad;
	return (glsl_noise::psrdnoise(p, {0.0f, 0.0f}, 0.0f, grad) + 1.0f) * 0.5f;
}

float SimplexNoise::Noise(glm::vec3 p)
{
	glm::vec3 grad;
	return (glsl_noise::psrdnoise(p, glm::vec3(0.0f, 0.0f, 0.0f), 0.0f, grad) +
			1.0f) *
		   0.5f;
}

float SimplexNoise::Noise(glm::vec4 p)
{
	return (glsl_noise::snoise(p) + 1.0f) * 0.5f;
}

float SimplexNoise::Noise2(glm::vec2 p)
{
	return CubicFactor(
		CubicFactor((Noise(p) + Noise(-p + glm::vec2{13, -27})) * 0.5f));
}

float SimplexNoise::Noise2(glm::vec3 p)
{
	return (Noise(p) + Noise(-p + glm::vec3{13, -27, 43})) * 0.5f;
}

float SimplexNoise::Noise2(glm::vec4 p)
{
	return (Noise(p) + Noise(-p + glm::vec4{13, -27, 43, -51})) * 0.5f;
}

float SimplexNoise::Fbm(glm::vec2 p, int octaves, float attenuation,
						float coordMultiplier, bool useGrad, bool useNoise2)
{
	p *= 5.0f;
	// 	float (SimplexNoise::*noise)(glm::vec2) = &SimplexNoise::Noise;
	// 	if (useNoise2) {
	// 		noise = &SimplexNoise::Noise2;
	// 	}

	glm::vec2 g{0, 0};
	// 	float dx = 0.001;

	float h = 0;
	float a = 1;
	float sum = 0;
	for (int i = 0; i < octaves; ++i) {
		float h0; //, h1, h2;

		glm::vec2 grad{0, 0};
		// 		h0 = (glsl_noise::perlin_cnoise(p) * 0.5 + 0.5) * a;
		h0 = (glsl_noise::psrdnoise(seedrot[i] * p, {0.0f, 0.0f}, i, grad) *
				  0.5 +
			  0.5) *
			 a;
		/*
		h0 = (this->*noise)(p)*a;
		*/
		if (useGrad) {
			/*
			h1 = (this->*noise)(p + glm::vec2{dx, 0.0f}) * a;
			h2 = (this->*noise)(p + glm::vec2{0.0f, dx}) * a;
			g += (glm::vec2{h1, h2} - h0) / dx;
			*/
			g += grad * a;
			h0 /= (1.0f + glm::dot(g, g));
		}
		h += h0;

		sum += a;
		p *= coordMultiplier;
		a *= attenuation;
	}
	return h / sum;
}

float SimplexNoise::NoiseRidges(glm::vec2 st)
{
	const float nscale = 5.0;
	glm::vec2 v = nscale * (st - 0.5f);
	const glm::vec2 p = glm::vec2(4.0, 4.0);
	float alpha = 0;
	glm::vec2 g, gsum;
	float warp = 0.13; // Nice "puffy clouds" warping

	float n = 0.0;
	float w = 1.0;
	float s = 1.0;
	gsum = glm::vec2(0.0);
	float sum = 0.0;
	for (int i = 0; i < 12; i++) {
		n += w * (glsl_noise::psrdnoise(seedrot[i + 17] * (s * v + warp * gsum),
										s * p, s * alpha, g) *
					  0.5 +
				  0.5);
		gsum += w * g;
		w *= 0.5;
		s *= 2.0;
		sum += w;
	}

	return 1.0f - n / sum;
}

float RecCell(glm::vec2 p, int octaves, float attenuation,
			  float coordMultiplier)
{
	p *= 5.0f;
	float h = 0;
	float a = 1;
	float sum = 0;
	for (int i = 0; i < octaves; ++i) {
		float h0;
		h0 = (glsl_noise::cellular(p).x * 0.5 + 0.5) * a;
		h += h0;

		sum += a;
		p *= coordMultiplier;
		a *= attenuation;
	}
	return h / sum;
}

float SimplexNoise::Terrain(glm::vec2 p, float verticalScale)
{
	// 	p /= 2.0f;
	// 	return RecCell(p, 14, 0.53, 1.7);
	// 	return glsl_noise::cellular(p).x;
	// 	glm::vec2 grad;
	// 	return NoiseRidges(p / 15.0f);
	// 	return glsl_noise::psrdnoise(p * 30.0f, {0.0f, 0.0f}, 0.0f, grad) * 0.5f
	// + 0.5f; 	return Fbm(p, 8, 3.0/7.0, 2.3, false, false) * (0.4 + 0.6 *
	// Noise(p * 0.1f + 123.f));
	int octaves = 14;
	bool useTwo = false;
	float biome = Noise(seedrot[33] * -p * 0.2f * 2.0f - 321.f);
	float scale = Noise(seedrot[34] * p * 0.4f + 123.f);
	float mountains = (NoiseRidges(p / 2.0f) + RecCell(p, 14, 0.53, 1.7) * 3) *
					  0.25 * (scale * 0.6 + 0.4);
	// 		Fbm(p, octaves, 0.43, 2.3, false, useTwo) * (scale * 0.6 + 0.4);
	float plains = Fbm(p * 0.3f - 100.0f, octaves, 0.5, 2.3, false, useTwo) *
				   (scale * 0.8 + 0.2);
	float h;
	if (biome < 0.25) {
		h = plains;
	} else if (biome < 0.5) {
		float f = (biome - 0.25) / 0.25;
		f = f * f * f * (f * (f * 6.0 - 15.0) + 10.0);
		h = plains + f * mountains;
	} else {
		h = plains + mountains;
	}
	return h * 0.5f;
}

} // namespace wg
