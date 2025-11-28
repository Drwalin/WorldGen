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
	bool parallel = false;

	constexpr static int OFF = 15;

	int iter = 0;
	
	void Init(int width, int height);
	Grid();
	~Grid();

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
};

#endif
