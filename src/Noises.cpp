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

#include <glm/matrix.hpp>

#include "../include/worldgen/Hash.hpp"
#include "../include/worldgen/Noises.hpp"

namespace wg
{

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

float textureBicubic(glm::vec2 p)
{

	// 	glm::vec2 texSize = textureSize(sampler, 0);
	// 	glm::vec2 invTexSize = 1.0 / texSize;

// 	p -= 0.5f;
// 
	glm::vec2 fxy = glm::fract(p);
	p -= fxy;

	glm::vec4 xc = cubic(fxy.x);
	glm::vec4 yc = cubic(fxy.y);

	glm::vec4 c =
		glm::vec4(p.x, p.x, p.y, p.y) + glm::vec4(-0.5, 1.5, -0.5, 1.5);

	glm::vec4 s = glm::vec4(glm::vec2(xc.x, xc.z) + glm::vec2(xc.y, xc.w),
							glm::vec2(yc.x, yc.z) + glm::vec2(yc.y, yc.w));
	glm::vec4 offset =
		c + glm::vec4(glm::vec2(xc.y, xc.w), glm::vec2(yc.y, yc.w)) / s;

	// 	offset *= invTexSize.xxyy;

	float sample0 = hx::Hash({offset.x, offset.z});
	float sample1 = hx::Hash({offset.y, offset.z});
	float sample2 = hx::Hash({offset.x, offset.w});
	float sample3 = hx::Hash({offset.y, offset.w});

	float sx = s.x / (s.x + s.y);
	float sy = s.z / (s.z + s.w);

	return glm::mix(glm::mix(sample3, sample2, sx),
					glm::mix(sample1, sample0, sx), sy);
}

static float sinc(float x) {
    const float PI = 3.1415926;
    return sin(PI * x) / (PI * x);
}

static float LWeight(float x) {
    if (fabs(x) < 1e-8) {
        return 1.0f;
    }
    else {
        return sinc(x) * sinc(x / 3);
    }
}

float Lanczos(glm::vec2 texCoord, int seed)
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
			float vecData =
				hx::Hash(samplePos+glm::vec2(m, n), seed);
			fY = LWeight(float(n) - interpFactor.y);
			float vecCoeef2 = fY;
			nSum += vecData * vecCoeef2 * vecCooef1;
		}
	}
	return nSum;
}

float Cubic(float a, float b, float f)
{
	float f2 = f*f;
	f = -2 * f * f2 + 3 * f2;
	return a * (1 - f) + b * f;
}

float Noise::NoiseV(glm::vec2 uv, int seed)
{
// 	return Lanczos(uv, seed);
// 	return textureBicubic(uv);

	glm::vec2 u = glm::floor(uv);
	glm::vec2 f = glm::fract(uv);
	float a = hx::Hash(u);
	float b = hx::Hash(u + glm::vec2(1.0, 0.0));
	float c = hx::Hash(u + glm::vec2(0.0, 1.0));
	float d = hx::Hash(u + glm::vec2(1.0, 1.0));
	return Cubic(Cubic(a, b, f.x), Cubic(c, d, f.x), f.y);
}

// 0: cubic
// 1: quintic
#define INTERPOLANT 0

// return value noise (in x) and its derivatives (in yzw)
float Noise::NoiseDV(glm::vec3 x)
{
	glm::ivec3 i = glm::ivec3(floor(x));
	glm::vec3 w = glm::fract(x);

#if INTERPOLANT == 1
	// quintic interpolation
	glm::vec3 u = w * w * w * (w * (w * 6.0f - 15.0f) + 10.0f);
	glm::vec3 du = 30.0f * w * w * (w * (w - 2.0f) + 1.0f);
#else
	// cubic interpolation
	glm::vec3 u = w * w * (3.0f - 2.0f * w);
#endif

	float a = hx::Hash(i + glm::ivec3(0, 0, 0));
	float b = hx::Hash(i + glm::ivec3(1, 0, 0));
	float c = hx::Hash(i + glm::ivec3(0, 1, 0));
	float d = hx::Hash(i + glm::ivec3(1, 1, 0));
	float e = hx::Hash(i + glm::ivec3(0, 0, 1));
	float f = hx::Hash(i + glm::ivec3(1, 0, 1));
	float g = hx::Hash(i + glm::ivec3(0, 1, 1));
	float h = hx::Hash(i + glm::ivec3(1, 1, 1));

	float k0 = a;
	float k1 = b - a;
	float k2 = c - a;
	float k3 = e - a;
	float k4 = a - b - c + d;
	float k5 = a - c - e + g;
	float k6 = a - b - e + f;
	float k7 = -a + b + c - d + e - f - g + h;

	return k0 + k1 * u.x + k2 * u.y + k3 * u.z + k4 * u.x * u.y +
		   k5 * u.y * u.z + k6 * u.z * u.x + k7 * u.x * u.y * u.z;
}

// return value noise (in x) and its derivatives (in yzw)
glm::vec4 Noise::NoiseD(glm::vec3 x)
{
	glm::ivec3 i = glm::ivec3(floor(x));
	glm::vec3 w = glm::fract(x);

#if INTERPOLANT == 1
	// quintic interpolation
	glm::vec3 u = w * w * w * (w * (w * 6.0f - 15.0f) + 10.0f);
	glm::vec3 du = 30.0f * w * w * (w * (w - 2.0f) + 1.0f);
#else
	// cubic interpolation
	glm::vec3 u = w * w * (3.0f - 2.0f * w);
	glm::vec3 du = 6.0f * w * (1.0f - w);
#endif

	float a = hx::Hash(i + glm::ivec3(0, 0, 0));
	float b = hx::Hash(i + glm::ivec3(1, 0, 0));
	float c = hx::Hash(i + glm::ivec3(0, 1, 0));
	float d = hx::Hash(i + glm::ivec3(1, 1, 0));
	float e = hx::Hash(i + glm::ivec3(0, 0, 1));
	float f = hx::Hash(i + glm::ivec3(1, 0, 1));
	float g = hx::Hash(i + glm::ivec3(0, 1, 1));
	float h = hx::Hash(i + glm::ivec3(1, 1, 1));

	float k0 = a;
	float k1 = b - a;
	float k2 = c - a;
	float k3 = e - a;
	float k4 = a - b - c + d;
	float k5 = a - c - e + g;
	float k6 = a - b - e + f;
	float k7 = -a + b + c - d + e - f - g + h;

	return glm::vec4(k0 + k1 * u.x + k2 * u.y + k3 * u.z + k4 * u.x * u.y +
						 k5 * u.y * u.z + k6 * u.z * u.x + k7 * u.x * u.y * u.z,
					 du * glm::vec3(k1 + k4 * u.y + k6 * u.z + k7 * u.y * u.z,
									k2 + k5 * u.z + k4 * u.x + k7 * u.z * u.x,
									k3 + k6 * u.x + k5 * u.y + k7 * u.x * u.y));
}

glm::vec4 Noise::FractalBrownianMotion(glm::vec3 x, int octaves)
{
	float f = 1.98; // could be 2.0f
	float s = 0.49; // could be 0.5
	float a = 0.0f;
	float b = 0.5;
	glm::vec3 d = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::mat3 m =
		glm::mat3(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	for (int i = 0; i < octaves; i++) {
		glm::vec4 n = NoiseD(x);
		a += b * n.x;						   // accumulate values
		d += b * m * glm::vec3(n.y, n.z, n.w); // accumulate derivatives
		b *= s;
		x = f * m /*3*/ * x;
		m = f * m /*3i*/ * m;
	}
	return glm::vec4(a, d);
}

float Noise::fbm(glm::vec3 x, float H, int octaves)
{
	float G = exp2(-H);
	float f = 1.0;
	float a = 1.0;
	float t = 0.0;
	for (int i = 0; i < octaves; i++) {
		t += a * Noise::NoiseD(f * x).x;
		f *= 2.0;
		a *= G;
	}
	return t;
}

float Noise::Terrain(glm::vec3 _p)
{
	constexpr glm::mat2 m = glm::mat2(0.8, -0.6, 0.6, 0.8);
	glm::vec2 p = {_p.x, _p.z};
	float a = 0.0f;
	float b = 1.0f;
	glm::vec2 d = glm::vec2(0.0f, 0.0f);
	for (int i = 0; i < 15; i++) {
		glm::vec3 n = NoiseD(glm::vec3(p.x, 0, p.y));
		d += glm::vec2(n.y, n.z);
		a += b * n.x / (1.0f + dot(d, d));
		b *= 0.5f;
		p = m * p * 2.0f;
	}
	return a;
}

float Noise::NoiseDV(glm::vec2 x)
{
	glm::vec2 i = floor(x);
	glm::vec2 f = fract(x);

	glm::vec2 u = f * f * f * (f * (f * 6.0f - 15.0f) + 10.0f);

	glm::vec2 ga = hx::Hash2(glm::vec3(i + glm::vec2(0.0f, 0.0f), 0.0f));
	glm::vec2 gb = hx::Hash2(glm::vec3(i + glm::vec2(1.0f, 0.0f), 0.0f));
	glm::vec2 gc = hx::Hash2(glm::vec3(i + glm::vec2(0.0f, 1.0f), 0.0f));
	glm::vec2 gd = hx::Hash2(glm::vec3(i + glm::vec2(1.0f, 1.0f), 0.0f));

	float va = dot(ga, f - glm::vec2(0.0f, 0.0f));
	float vb = dot(gb, f - glm::vec2(1.0f, 0.0f));
	float vc = dot(gc, f - glm::vec2(0.0f, 1.0f));
	float vd = dot(gd, f - glm::vec2(1.0f, 1.0f));

	return va + u.x * (vb - va) + u.y * (vc - va) +
		   u.x * u.y * (va - vb - vc + vd);
}

// return value noise (in x) and its derivatives (in yzw)
glm::vec3 Noise::NoiseD(glm::vec2 x)
{
	glm::vec2 i = floor(x);
	glm::vec2 f = fract(x);

	glm::vec2 u = f * f * f * (f * (f * 6.0f - 15.0f) + 10.0f);
	glm::vec2 du = 30.0f * f * f * (f * (f - 2.0f) + 1.0f);

	glm::vec2 ga = hx::Hash2(glm::vec3(i + glm::vec2(0.0f, 0.0f), 0.0f));
	glm::vec2 gb = hx::Hash2(glm::vec3(i + glm::vec2(1.0f, 0.0f), 0.0f));
	glm::vec2 gc = hx::Hash2(glm::vec3(i + glm::vec2(0.0f, 1.0f), 0.0f));
	glm::vec2 gd = hx::Hash2(glm::vec3(i + glm::vec2(1.0f, 1.0f), 0.0f));

	float va = dot(ga, f - glm::vec2(0.0f, 0.0f));
	float vb = dot(gb, f - glm::vec2(1.0f, 0.0f));
	float vc = dot(gc, f - glm::vec2(0.0f, 1.0f));
	float vd = dot(gd, f - glm::vec2(1.0f, 1.0f));

	return glm::vec3(va + u.x * (vb - va) + u.y * (vc - va) +
						 u.x * u.y * (va - vb - vc + vd), // value
					 ga + u.x * (gb - ga) + u.y * (gc - ga) +
						 u.x * u.y * (ga - gb - gc + gd) + // derivatives
						 du * ((glm::vec2(u.y, u.x)) * (va - vb - vc + vd) +
							   glm::vec2(vb, vc) - va));
}

glm::vec3 Noise::FractalBrownianMotion(glm::vec2 x, int octaves)
{
	float f = 1.98; // could be 2.0f
	float s = 0.49; // could be 0.5
	float a = 0.0f;
	float b = 0.5;
	glm::vec3 d = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::mat3 m =
		glm::mat3(1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	for (int i = 0; i < octaves; i++) {
		float h0 = NoiseV(x);
		float h1 = NoiseV(x + glm::vec2(0.0625f, 0));
		float h2 = NoiseV(x + glm::vec2(0, 0.0625f));
		glm::vec3 n {
			h0,
			-(h1 - h0) / 0.0625f,
			-(h2 - h0) / 0.0625f};
		
		
		
		a += b * n.x;					 // accumulate values
		d += b * glm::vec3(n.y, n.z, 0); // accumulate derivatives
		b *= s;
		x = f * m /*3*/ * glm::vec3(x, 0);
		m = f * m /*3i*/ * m;
	}
	return glm::vec3(a, d.x, d.y);
}

float Noise::fbm(glm::vec2 x, float H, int octaves)
{
	return FractalBrownianMotion(x, octaves).x;
	
	float G = exp2(-H);
	float f = 1.0;
	float a = 1.0;
	float t = 0.0;
	for (int i = 0; i < octaves; i++) {
		t += a * Noise::NoiseV(f * x);
		f *= 2.0;
		a *= G;
	}
	return t;
}


float Noise::Ridges(glm::vec2 x, float horizontalScale, float steep, float octaves)
{
	float h = 0;
	float f = 1;
	glm::vec2 grad = {0,0};
	
	float div = 1.0f;
	
	for (int i = 0; i < octaves; i++) {
		float h0 = NoiseV(x, i);
		float h1 = NoiseV(x + glm::vec2(0.0625f, 0), i);
		float h2 = NoiseV(x + glm::vec2(0, 0.0625f), i);
		grad += glm::vec2{
			(h1 - h0) / (0.0625f*horizontalScale),
			(h2 - h0) / (0.0625f*horizontalScale)};
		
		float gf = glm::dot(grad, grad) * 0.7f;
		gf = 1.0f / (1.0f + gf);
		h += h0 * gf * f;
		f *= steep;
		
		div += f;
		x *= 1.3f;
	}
	return h / div;
}

float Noise::Terrain(glm::vec2 p, float horizontalScale)
{
	return Ridges(p, horizontalScale, 0.85, 15);
}


} // namespace wg
