
#pragma once

#ifndef HYDRO_EROSION_HPP
#define HYDRO_EROSION_HPP

#include <cstdio>
#include <cmath>
#include <cinttypes>

struct Flux {
	float L;
	float T;
	float R;
	float B;
};

struct FloatArray4 {
	float v[4];
};

using TileId = int32_t;

class HydroErosion {
public:
	
	union {
		float *b;
		float *ground;
	};
	union {
		float *water;
		float *d;
	};
	union {
		float *sediment;
		float *suspendedSediment;
		float *s;
	};
	union {
		float *newSediment;
	};
	union {
		Flux *f;
		FloatArray4 *fluxArray;
	};
	float *vx;
	float *vy;
	
public:
	
	bool wrap;
	
	const int width, height;
	float dt;
	union {
		float softness;
		float Ks;
		float dissolvingConstant;
	};
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
	
	HydroErosion(int width, int height, float l) : width(width), height(height), l(l) {
		wrap = false;//true;
		dt = 0.02;
		A = l*l;
		g = 9.81;
		Kd = 0.1;
		Kc = 0.5;
		minimumSedimentCapacity = 0.1;
		
		ground = new float[width*height];
		water = new float[width*height];
		suspendedSediment = new float[width*height];
		newSediment = new float[width*height];
		Ks = 0.1;
		f = new Flux[width*height];
		vx = new float[width*height];
		vy = new float[width*height];
		for(int i=0; i<width*height; ++i) {
			ground[i] = 0;
			water[i] = 0;
			suspendedSediment[i] = 0;
			vx[i] = 0;
			vy[i] = 0;
		}
	}
	
	~HydroErosion() {
		delete ground;
		delete water;
		delete suspendedSediment;
		delete newSediment;
		delete f;
		delete vx;
		delete vy;
	}
	
public:
	
	template<bool safe>
	TileId At(int x, int y) const;
	template<bool safe, int dir>
	TileId Neighbour(int x, int y) const;
	template<bool safe>
	TileId Neighbour(int x, int y, int dir) const;
	
	template<bool safe, int dir>
	float CalcFluxInDirection(TileId src, TileId neigh) const;
	void LimitFlux(TileId src);
	template<bool safe>
	void CalcOutflux(int x, int y); // 3.2.1
	
	template<bool safe>
	void UpdateWaterLevelAndVelocity(int x, int y); // 3.2.2

	
	template<bool safe>
	float SinusLocalTiltAngle(TileId t, int x, int y) const;
	template<bool safe>
	void ErosionAndDeposition(int x, int y); // 3.3
	
	template<bool safe>
	void SedimentTransportation(int x, int y); // 3.4
	template<bool safe>
	void SedimentTransportationUpdateWhatsLeft(int x, int y); // 3.4
	
	float EvaporationRate(int x, int y);
	template<bool safe>
	void Evaporation(int x, int y); // 3.5
	
	// to be executed after water increase
	void FullCycle();
};

struct Tile {
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
};

struct Grid {
	Grid() {
		width = height = 0;
		tiles = NULL;
		dt = 0.01;
		crossSectionalAreaOfPipe = 1;
		gravity = 9.81;
		tileDimensionSize = 1;
		
		depositionConstant = 1;//0.01;
		sedimentCapacityConstant = 1;//0.5;
		minimumSedimentCapacity = 0.1;
	}
	~Grid() {
		if(tiles)
			delete[] tiles;
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
	inline float CalcFluxInDirection(Tile& src, Tile& neigh) const;
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
	inline void ErosionAndDeposition(int x, int y); // 3.3
	
	template<bool safe>
	inline void SedimentTransportation(int x, int y); // 3.4
	
	inline float EvaporationRate(int x, int y);
	template<bool safe>
	inline void Evaporation(int x, int y); // 3.4
	
	// to be executed after water increase
	inline void FullCycle();
};

#include "HydroErosion.incl.hpp"

#endif

