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

struct Flux {
	float L;
	float B;
	float R;
	float T;
};

struct FluxUnion {
	FluxUnion() {
		f.L = f.B = f.R = f.T = 0;
	}
	union {
		Flux f;
		float fluxArray[4];
	};
};

struct Velocity {
	float x = 0.0f, y = 0.0f;
};

struct Grid {
	int iter = 0;
	void Init(int width, int height) {
		this->width = width;
		this->height = height;
		ground = new float[width*height] - 1;
		water = new float[width*height] - 1;
		sediment = new float[width*height] - 1;
		hardness = new float[width*height] - 1;
		deltaSedimentGround = new float[width*height] - 1;
		for (int i=1; i<=width*height; ++i) {
			ground[i] = 0.0f;
			water[i] = 0.0f;
			sediment[i] = 0.0f;
			hardness[i] = 0.1f;
			deltaSedimentGround[i] = 0.0f;
		}
		velocity = new Velocity[width*height] - 1;
		flux = new FluxUnion[width*height] - 1;
		
	}
	Grid() {
		width = height = 0;
		dt = 0.02;
		crossSectionalAreaOfPipe = .6;
		gravity = 9.81;
		tileDimensionSize = 1;
		
		depositionConstant = 0.03;
		sedimentCapacityConstant = 0.03;
		minimumSedimentCapacity = 0.1 * 0;
	}
	~Grid() {
		if (ground) { delete[] ground; }
		if (water) { delete[] water; }
		if (sediment) { delete[] sediment; }
		if (hardness) { delete[] hardness; }
		if (deltaSedimentGround) { delete[] deltaSedimentGround; }
		if (velocity) { delete[] velocity; }
		if (flux) { delete[] flux; }
	}
	
	union {
		float *b = nullptr;
		float *ground;
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
		float *hardness = nullptr;
		float *Ks;
		float *dissolvingConstant;
	};
	float *deltaSedimentGround = nullptr;
	Velocity *velocity = nullptr;
	FluxUnion *flux = nullptr;
	
	
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

