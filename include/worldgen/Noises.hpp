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

#include "../../thirdparty/OpenSimplex2/OpenSimplex2F/OpenSimplex2F.hpp"

namespace wg
{
class Noise
{
public:
	static float NoiseV(glm::vec2 x, int seed = 0);
	static glm::vec3 NoiseVG(glm::vec2 x, int seed = 0);
	
	static float NoiseDV(glm::vec3 x);
	static float NoiseDV(glm::vec2 x);
	
	static glm::vec4 NoiseD(glm::vec3 x);
	static glm::vec4 FractalBrownianMotion(glm::vec3 x, int octaves);
	static float Terrain(glm::vec3 p);
	static float fbm(glm::vec3 x, float H, int octaves);
	
	static glm::vec3 NoiseD(glm::vec2 x);
	static glm::vec3 FractalBrownianMotion(glm::vec2 x, int octaves);
	static float fbm(glm::vec2 x, float H, int octaves);
	
	static float Ridges(glm::vec2 x, float horizontalScale, float steep, float octaves, float gradientInfluence, int seed);
	static float Ridges2(glm::vec2 x, float horizontalScale, float steep, float octaves, float gradientInfluence, int seed);
	
	static float Terrain(glm::vec2 p, float horizontalScale);
};

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
	
	float Fbm(glm::vec2 p, int octaves, float attenuation, float coordMultiplier, bool useGrad, bool useNoise2, float verticalScale);
	
	float Terrain(glm::vec2 p, float verticalScale);
	
private:
	static OpenSimplex2F::OpenSimplexEnv *ose;
	OpenSimplex2F::OpenSimplexGradients *osg = nullptr;
};
} // namespace wg
