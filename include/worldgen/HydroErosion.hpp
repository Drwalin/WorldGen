#pragma once

#ifndef HYDRO_EROSION_HPP
#define HYDRO_EROSION_HPP

#include <cstdio>

struct Flux {
	float L;
	float B;
	float R;
	float T;
};

struct alignas(64) Tile {
	Tile() {
		ground = 0;
		water = 0;
		sediment = 0;
		f.L = f.B = f.R = f.T = 0;
		vx = vy = 0;
		
		hardness = 0.1;
	}
	union {
		float b;
		float ground;
	};
	union {
		float water;
		float d;
	};
	union {
		float sediment;
		float suspendedSediment;
		float s;
	};
	union {
		float hardness;
		float Ks;
		float dissolvingConstant;
	};
	union {
		Flux f;
		float fluxArray[4];
	};
	float vx, vy;
	float deltaSedimentGround;
};

struct Grid {
	int iter = 0;
	void Init(int width, int height) {
		this->width = width;
		this->height = height;
		this->tiles = new Tile[width*height];
	}
	Grid() {
		width = height = 0;
		tiles = NULL;
		dt = 0.1;
		crossSectionalAreaOfPipe = .6;
		gravity = 9.81;
		tileDimensionSize = 1;
		
		depositionConstant = 0.03;
		sedimentCapacityConstant = 0.03;
		minimumSedimentCapacity = 0.1 * 0;
	}
	~Grid() {
		if(tiles) {
			delete[] tiles;
		}
	}
	Tile* tiles;
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
	inline Tile* At(int x, int y) const;
	template<bool safe, int dir>
	inline Tile* Neighbour(int x, int y) const;
	template<bool safe>
	inline Tile* Neighbour(int x, int y, int dir) const;
	
	template<bool safe, int dir>
	inline float CalcFluxInDirection(Tile& src, Tile* neigh) const;
	inline void LimitFlux(Tile& src);
	template<bool safe>
	void CalcOutflux(int x, int y); // 3.2.1
	
	template<bool safe>
	void UpdateWaterLevel(Tile& src, Tile** neighs);
	template<bool safe>
	void UpdateWaterLevelAndVelocity(int x, int y); // 3.2.2
	
	template<bool safe>
	inline float SinusLocalTiltAngle(Tile& t, int x, int y);
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

	template<int BORDER, bool PARALLEL, typename T1, typename T2>
	inline void ForEachSafeBorders(T1 &&funcSafe, T2 &&funcUnsafe);
	
	// to be executed after water increase
	void FullCycle();
};

#include "HydroErosion.incl.hpp"
#ifdef HYDRO_EROSION_INCL_HPP
#endif

#endif

