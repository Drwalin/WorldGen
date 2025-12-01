#pragma once

#ifndef HYDRO_EROSION_HPP
#define HYDRO_EROSION_HPP

#include "../../../OpenGLWrapper/include/openglwrapper/VBO.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Shader.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Texture.hpp"

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
	bool useGpu;
	bool useWater = true;
	bool useThermalErosion = true;
	bool useSmoothing = false;
	bool parallel = false;

	constexpr static int OFF = 15;

	int iteration = 0;
	
	void Init(int width, int height, bool useGpu);
	Grid();
	~Grid();

	GroundLayers *ground = nullptr;
	float *water = nullptr;
	float *sediment = nullptr;
	float *temp1 = nullptr;
	Velocity *velocity = nullptr;
	Flux *flux = nullptr;

	float hardness[2] = {0.003, 0.02};
// 	float hardness[2] = {0.01, 0.04};
	// tan(30) ~= 0.577
	// tan(45) ~= 1
	// tan(60) ~= 1.732
	// tan(75) ~= 3.732
	static constexpr float tangentOfAngleOfRecluse[2] = {1.7f, 1.0f};
	// {(float)tan(M_PI/4.0f), (float)tan(M_PI/6.0f)}; // 60*, 45*

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

	void CallHydroErosion();
	void CallThermalErosion();
	void CallSmoothing();
	void CallEvaporation();

	template <typename TFunc> void ForEachSafeBorders(TFunc &&funcSafe);

	float SumFlux(int t);

	// to be executed after water increase
	void FullCycle();
	void UpdateHeightsTexture(gl::Texture *tex);
	
	struct GPUCompute {
		~GPUCompute();
		
		gl::Shader shaderCalcOutFlux;
		gl::Shader shaderUpdateWaterLevelAndVelocity;
		gl::Shader shaderErosionAndDepositionCalculation;
		gl::Shader shaderErosionAndDepositionUpdate;
		gl::Shader shaderSedimentTransportation;
		gl::Shader shaderSedimentTransportationUpdate;
		
		gl::Shader shaderThermalErosionCalculation;
		gl::Shader shaderThermalErosionUpdate;
		
		gl::Shader shaderEvaporation;
		
		gl::Shader shaderSmooth;
		gl::Shader shaderSmoothUpdate;
		
		gl::Shader shaderUpdateRainAndRiver;
		
		gl::Shader shaderUpdateHeightTexture;
		int textureUniformLocation;
		
		void CallCalcOutFlux();
		void CallUpdateWaterLevelAndVelocity();
		void CallErosionAndDepositionCalculation();
		void CallErosionAndDepositionUpdate();
		void CallSedimentTransportation();
		void CallSedimentTransportationUpdate();
		void CallThermalErosionCalculation();
		void CallThermalErosionUpdate();
		void CallEvaporation();
		void CallSmooth();
		void CallSmoothUpdate();
		
		void CallShader(gl::Shader *shader);
		
		void CallUpdateRainAndRiver();
		void UpdateHeightsTexture(gl::Texture *tex);
		
		
		gl::VBO *vboGround = nullptr;
		gl::VBO *vboWater = nullptr;
		gl::VBO *vboSediment = nullptr;
		gl::VBO *vboTemp1 = nullptr;
		gl::VBO *vboVelocity = nullptr;
		gl::VBO *vboFlux = nullptr;
		
		gl::VBO *vboRiverSources = nullptr;
		
		int width, height;
		
		void UpdateGround(GroundLayers *data);
		void UpdateWater(float *data);
		void UpdateSediment(float *data);
		void UpdateTemp1(float *data);
		void UpdateVelocity(Velocity *data);
		void UpdateFlux(Flux *data);
		void UpdateRiverSources(glm::vec3 *data, int amount);
		
		
		void BindBuffers();
		void SetUniforms(gl::Shader *shader);
		
		void Init(int w, int h, Grid *grid);
		
		Grid *grid;
	} gpu;
};

#endif
