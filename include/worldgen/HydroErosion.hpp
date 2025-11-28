#pragma once

#ifndef HYDRO_EROSION_HPP
#define HYDRO_EROSION_HPP

/*
 * Angles of repose:
 *    dry sand: 34*
 *    wet sand: 45*
 *    soil: 30* - 45*
 *    snow: 38*
 *    sand with water: 15* - 30*    // maybe can be used for sediment
 */

#if not defined EROSION_STRUCTS_DEFINED
#define EROSION_STRUCTS_DEFINED
struct Flux {
	float f[4];
};

struct Velocity {
	float x;
	float y;
};

struct GroundLayers {
	float layers[2];
};
#endif

struct Grid {
	bool useWater = true;
	bool useThermalErosion = true;
	bool useSmoothing = false;
	
	constexpr static int OFF = 15;
	
	int iter = 0;
	void Init(int width, int height) {
		this->width = width;
		this->height = height;
		ground = new GroundLayers[width*height + OFF + 1] + OFF;
		water = new float[width*height + OFF + 1] + OFF;
		sediment = new float[width*height + OFF + 1] + OFF;
		deltaSedimentGround = new float[width*height + OFF + 1] + OFF;
		velocity = new Velocity[width*height + OFF + 1] + OFF;
		flux = new Flux[width*height + OFF + 1] + OFF;
		for (int i=1; i<=width*height; ++i) {
			water[i] = 0.0f;
			sediment[i] = 0.0f;
			deltaSedimentGround[i] = 0.0f;
			flux[i].f[0] = 0;
			flux[i].f[1] = 0;
			flux[i].f[2] = 0;
			flux[i].f[3] = 0;
			ground[i].layers[0] = 0;
			ground[i].layers[1] = 0;
		}
	}
	Grid() {
		width = height = 0;
		dt = 0.03;
		crossSectionalAreaOfPipe = .6;
		gravity = 9.81;
		tileDimensionSize = 1;
		
		depositionConstant = 0.03;
		sedimentCapacityConstant = 0.03;
		minimumSedimentCapacity = 0.1;
	}
	~Grid() {
		if (ground) { delete[] (ground - OFF); }
		if (water) { delete[] (water - OFF); }
		if (sediment) { delete[] (sediment - OFF); }
		if (deltaSedimentGround) { delete[] (deltaSedimentGround - OFF); }
		if (velocity) { delete[] (velocity - OFF); }
		if (flux) { delete[] (flux - OFF); }
	}
	
	union {
		GroundLayers *b = nullptr;
		GroundLayers *ground;
	};
	union {
		float *water = nullptr;
		float *d;
	};
	union {
		float *sediment = nullptr;
		float *suspendedSediment;
		float *s;
	};
	union {
		float hardness[2] = {0.03, 0.07};
		float Ks[2];
		float dissolvingConstant[2];
	};
	float *deltaSedimentGround = nullptr;
	Velocity *velocity = nullptr;
	Flux *flux = nullptr;
	
	
	// tan(30) ~= 0.577
	// tan(45) ~= 1
	// tan(60) ~= 1.732
	// tan(75) ~= 3.732
	static constexpr float tangentOfAngleOfRecluse[2] = {1.7f, 1.0f};
	// {(float)tan(M_PI/4.0f), (float)tan(M_PI/6.0f)}; // 60*, 45*
	
	
	int width, height;
	float dt;
	union {
		float A;
		float crossSectionalAreaOfPipe;
	};
	union {
		float g;
		float gravity;
	};
	union {
		float l;
		float tileDimensionSize;
	};
	union {
		float Kd;
		float depositionConstant;
	};
	union {
		float Kc;
		float sedimentCapacityConstant;
	};
	float minimumSedimentCapacity;
	
	template<bool safe>
	inline int At(int x, int y) const;
	template<bool safe, int dir>
	inline int Neighbour(int x, int y) const;
	template<bool safe>
	inline int Neighbour(int x, int y, int dir) const;
	
	template<bool safe, int dir>
	inline float CalcFluxInDirection(int src, int neigh) const;
	inline void LimitFlux(int src);
	template<bool safe>
	void CalcOutflux(int x, int y); // 3.2.1
	
	template<bool safe>
	void UpdateWaterLevel(int src, int* neighs);
	template<bool safe>
	void UpdateWaterLevelAndVelocity(int x, int y); // 3.2.2
	
	template<bool safe>
	inline float SinusLocalTiltAngle(int t, int x, int y);
	template<bool safe>
	inline void ErosionAndDepositionCalculation(int x, int y); // 3.3
	template<bool safe>
	inline void ErosionAndDepositionUpdate(int x, int y); // 3.3
	
	template<bool safe>
	inline void SedimentTransportation(int x, int y); // 3.4
	template<bool safe>
	inline void SedimentTransportationUpdate(int x, int y); // 3.4
	
	inline float EvaporationRate(int x, int y);
	template<bool safe>
	inline void Evaporation(int x, int y); // 3.4
	template<bool safe>
	inline void Smooth(int x, int y); // 3.4
	template<bool safe>
	inline void SmoothUpdate(int x, int y); // 3.4
	template<bool safe>
	inline void ClearDelta(int x, int y);
	
	
	template<bool safe>
	inline void ThermalErosionCalculation(int x, int y);
	template<bool safe>
	inline void ThermalErosionUpdate(int x, int y);

	template<int BORDER, bool PARALLEL, typename T1, typename T2>
	inline void ForEachSafeBorders(T1 &&funcSafe, T2 &&funcUnsafe);
	
	
	inline float SumFlux(int t);
	
	// to be executed after water increase
	void FullCycle();
};

#include "HydroErosion.incl.hpp"
#ifdef HYDRO_EROSION_INCL_HPP
#endif

#endif

