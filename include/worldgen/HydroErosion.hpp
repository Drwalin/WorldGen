#pragma once

#ifndef HYDRO_EROSION_HPP
#define HYDRO_EROSION_HPP

#include <vector>
#include <functional>

#include "../../../OpenGLWrapper/include/openglwrapper/VBO.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Shader.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Texture.hpp"

namespace gl
{
class VBO;
class Shader;
class TExture;
} // namespace gl

/*
 * Angles of repose:
 *    dry sand: 34*
 *    wet sand: 45*
 *    soil: 30* - 45*
 *    snow: 38*
 *    sand with water: 15* - 30*    // maybe can be used for sediment
 */

#include "HydroErosionStructs.h"

struct Grid {
	int parallelThreads = 1;
	bool useGpu;

	bool useWater = true;
	bool useThermalErosion = !true;
	bool useSmoothing = !true;

	bool parallel = false;

	constexpr static int ALIGNEMENT = 16;
	constexpr static int PADDING = ALIGNEMENT - 1;
	constexpr static int OFFSET = ALIGNEMENT - 1;

	int iteration = 0;

	void Init(int width, int height, bool useGpu, int workGroupSizeDim=16);
	Grid();
	~Grid();

	uint32_t elements;
	uint32_t elementsStorage;

	glm::vec4 *water_sediment_temp = nullptr;

	GroundLayers *ground = nullptr;
	Velocity *velocity = nullptr;
	Flux *flux = nullptr;
	float *water = nullptr;
	float *sediment = nullptr;
	float *temp1 = nullptr;
	float *temp2 = nullptr;

	// float hardness[2] = {0.008, 0.02};
	// float hardness[2] = {0.01, 0.04};
	// float hardness[2] = {0.04, 0.2};
	float hardness[2] = {0.1, 0.5};
	// float hardness[2] = {0.3, 0.8};
	// tan(30) ~= 0.577
	// tan(34) ~= 0.675
	// tan(45) ~= 1
	// tan(60) ~= 1.732
	// tan(75) ~= 3.732
	float tangentOfAngleOfRecluse[2] = {1, 0.675}; // 45* 34*
	// float tangentOfAngleOfRecluse[2] = {1.7f, 1.0f}; // 60* 45*

	int width, height;
	float dt = 0.03;
	union {
		float A = 0.6;
		float crossSectionalAreaOfPipe;
	};
	union {
		float g;
		float gravity = 9.81;
	};
	union {
		float l = 1.0;
		float tileDimensionSize;
	};
	union {
		float Kd = 0.03;
		float depositionConstant;
	};
	union {
		float Kc = 0.03;
		float sedimentCapacityConstant;
	};
	float minimumSedimentCapacity = 0.03;

	inline int At(int x, int y) const
	{
		if (x < 0)
			return 0;
		else if (x >= width)
			return 0;
		if (y < 0)
			return 0;
		else if (y >= height)
			return 0;
		return (x * height + y) + 1;
	}

	template <typename TFunc> void ForEachSafeBorders(TFunc &&funcSafe);

	float SumFlux(int t);

	// to be executed after water increase
	void FullCycle();
	void UpdateHeightsTexture(gl::Texture *tex);

	struct StageData {
		std::function<void()> functionCpu;
		std::function<void(gl::Shader *)> functionGpu;
		gl::Shader *shader = nullptr;
		std::string functionName;

		~StageData() { delete[] shader; }
	};
	std::vector<std::vector<StageData>> stages;

	struct GPUCompute {
		~GPUCompute();

		gl::Shader *shaderUpdateHeightTexture;
		int textureUniformLocation;

		void CallShader(gl::Shader *shader);

		void UpdateHeightsTexture(gl::Texture *tex);

		gl::VBO *vboGround = nullptr;
		gl::VBO *vboVelocity = nullptr;
		gl::VBO *vboFlux = nullptr;
		gl::VBO *vboWaterSedimentTemp = nullptr;

		int width, height;

		void UpdateAll();
		void UpdateGround(const GroundLayers *data);
		void UpdateWater(const float *data);
		void UpdateSediment(const float *data);
		void UpdateTemp1(const float *data);
		void UpdateTemp2(const float *data);
		void UpdateVelocity(const Velocity *data);
		void UpdateFlux(const Flux *data);

		void FetchAll();
		void FetchGround(GroundLayers *data);
		void FetchWater(float *data);
		void FetchSediment(float *data);
		void FetchTemp1(float *data);
		void FetchTemp2(float *data);
		void FetchVelocity(Velocity *data);
		void FetchFlux(Flux *data);

		void BindBuffers();
		void SetUniforms(gl::Shader *shader);

		void Init(int w, int h, Grid *grid, int workGroupSizeDim);

		Grid *grid;
	} gpu;
};

#endif
