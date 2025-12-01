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

#pragma once

#include <glm/common.hpp>
#include <glm/fwd.hpp>
#include <glm/vec2.hpp>
#include <glm/mat2x2.hpp>

namespace wg
{
class SimplexNoise
{
public:
	SimplexNoise(uint64_t seed=12345);
	~SimplexNoise();
	
	void Init(uint64_t seed);
	
	float Noise(glm::vec2 p);
	float Noise(glm::vec3 p);
	float Noise(glm::vec4 p);
	
	float Noise2(glm::vec2 p);
	float Noise2(glm::vec3 p);
	float Noise2(glm::vec4 p);
	
	float NoiseRidges(glm::vec2 st);
	float Fbm(glm::vec2 p, int octaves, float attenuation, float coordMultiplier, bool useGrad, bool useNoise2);
	
	float Terrain(glm::vec2 p, float verticalScale);
	
private:
	uint64_t seed;
	static constexpr int AMOUNT_SEEDED_VALS = 64;
	float seedf[AMOUNT_SEEDED_VALS];
	float seedrad[AMOUNT_SEEDED_VALS];
	uint64_t seedi[AMOUNT_SEEDED_VALS];
	glm::mat2 seedrot[AMOUNT_SEEDED_VALS];
};
} // namespace wg
