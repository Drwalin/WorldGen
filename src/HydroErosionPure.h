
#include "../include/worldgen/GlslSharedWithCpp_Header.h"
#include "../include/worldgen/HydroErosionMacros.h"
#include "../include/worldgen/HydroErosionStructs.h"

#if GL_core_profile
#else
using uint = unsigned int;
namespace glsl_noise
{
using namespace glm;
#endif
#include "../thirdparty/webgl-noise/src/cellular2D.glsl"
#include "../thirdparty/psrdnoise/src/psrdnoise3.glsl"
#include "../thirdparty/psrdnoise/src/psrdnoise2.glsl"
#if GL_core_profile
#else
}
#endif

#include "./CGlslRandom.h"

#if GL_core_profile
#else
namespace HydroPure
{
using namespace glm;
using namespace glsl_noise;
#endif

static uniform float hardness[2] = {0.03, 0.07};
static uniform float tangentOfAngleOfRecluse[2] = {1.7, 1};
static uniform int width, height;
static uniform float dt;
static uniform float A;	 // crossSectionalAreaOfPipe
static uniform float g;	 // gravity
static uniform float l;	 // tileDimensionSize
static uniform float Kd; // depositionConstant
static uniform float Kc; // sedimentCapacityConstant
static uniform float minimumSedimentCapacity;
static uniform int iteration;

layout(std430, binding = 1) BUFFER(BufferBlock1, GroundLayers, ground);
layout(std430, binding = 2) BUFFER(BufferBlock2, Velocity, velocity);
layout(std430, binding = 3) BUFFER(BufferBlock3, Flux, flux);

#if GL_core_profile
layout(std430, binding = 4) buffer BufferBlock7 {
	float paddingWater[PADDING];
	float water[ELEMENTS];
	
	float paddingSediment[PADDING];
	float sediment[ELEMENTS];
	
	float paddingTemp1[PADDING];
	float temp1[ELEMENTS];
	
// 	float paddingTemp2[PADDING];
// 	float temp2[ELEMENTS];
};
#else
static float *water = nullptr;
static float *sediment = nullptr;
static float *temp1 = nullptr;
#endif

#if GL_core_profile
#else
namespace _____holder
{
inline void _holder() { (void)iteration; }
} // namespace _____holder
#endif

inline float TotalGround(int t)
{
	GroundLayers g = ground[t];
	return g.layers[0] + g.layers[1];
}

inline GroundLayers AddGeneralGround(GroundLayers g, float dv)
{
	g.layers[1] += dv;
	if (g.layers[1] < 0) {
		g.layers[0] += g.layers[1];
		g.layers[1] = 0;
	}
	return g;
}

inline int At(int x, int y)
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

inline int Neighbour(int x, int y, int dir)
{
	switch (dir) {
	case 0:
		return At(x - 1, y);
	case 1:
		return At(x, y + 1);
	case 2:
		return At(x + 1, y);
	case 3:
		return At(x, y - 1);
	}
	return 0;
}

struct RiverSource {
	ivec2 coord;
	float amount;
};

inline RiverSource CalcRiverSource(ivec2 p)
{
	const uint sourcesGridSize = 97u;
	uvec2 v = uvec2(p);
	uvec2 md = v % sourcesGridSize;
	uvec2 base = v / sourcesGridSize;

	uvec4 rnd = RandomUint(ivec3(int(base.x), int(base.y), 0));

	uvec2 source = v - md + (uvec2(rnd.x, rnd.y) % sourcesGridSize);
	if (rnd.z % 7 == 0) {
		return RiverSource(ivec2(source), 0.0);
	}
	float amount = float(rnd.w % 1511u) / 1024.0 + 0.3;
	amount = amount * amount;
	return RiverSource(ivec2(source), amount);
}

inline void UpdateRainAndRiver(int x, int y)
{
	int src = At(x, y);
	float w = water[src];

	vec3 coord =
		vec3(float(x), float(y), float(iteration % 100000) * 0.02) * float(0.01);
	vec3 grad;
	float val = (psrdnoise(coord, vec3(0, 0, 0),
						  0.001 * float(iteration % 100000), grad) *
					0.5 +
				0.5);
	val = clamp(val, float(0.0), float(1.0));
	val *= val;
	val *= val;
	val *= val;
	val *= 0.01;

	RiverSource rs = CalcRiverSource(ivec2(x, y));
	ivec2 ap = (rs.coord - ivec2(x, y));
	int dp = ap.x * ap.x + ap.y * ap.y;
	if (dp <= 1) {
		val += rs.amount / (1.0 + dp);// * (1.0 / 5.0);
	}

	water[src] = w + val * dt;
}

inline float SumFluxF(Flux f) { return f.f[0] + f.f[1] + f.f[2] + f.f[3]; }

inline float SumFlux(int t) { return SumFluxF(flux[t]); }

inline float CalcFluxInDirection(int src, int neigh, int dir, Flux fl,
								 float srcSum)
{
	if (neigh == 0)
		return 0;
	int dst = neigh;
	float dh = srcSum - (TotalGround(dst) + sediment[dst] + water[dst]);
	float f = fl.f[dir] + dt * A * g * dh / l;
	if (f < 0)
		f = 0;
	return f;
}

inline Flux LimitFlux(int src, Flux f, float waterLevel)
{
	const float outflux = SumFluxF(f);
	const float _water = waterLevel * l * l;
	if (outflux <= 0.001) {
		f.f[0] = 0;
		f.f[1] = 0;
		f.f[2] = 0;
		f.f[3] = 0;
		return f;
	}
	float K = _water / (outflux * dt);
	if (K > 1) {
		K = 1;
	}
	if (K >= 0) {
		FOR_EACH_DIR({ f.f[DIR] *= K; })
	}
	return f;
}

inline void CalcOutFlux(int x, int y)
{
	int src = At(x, y);
	Flux f = flux[src];
	float waterLevel = water[src];
	float srcSum = TotalGround(src) + sediment[src] + waterLevel;
	NEIGHBOURS(neighs, x, y);
	FOR_EACH_DIR_COND(
		neighs[DIR] != 0,
		(f.f[DIR] = CalcFluxInDirection(src, neighs[DIR], DIR, f, srcSum));)
	flux[src] = LimitFlux(src, f, waterLevel);
}

inline float UpdateWaterLevel(Flux srcFlux, int neighs[4], float oldWater,
							  Flux neighFlux[4])
{
	float fs = 0;
	FOR_EACH_DIR_COND(neighs[DIR] != 0,
					  fs += neighFlux[DIR].f[R_DIR] - srcFlux.f[DIR];)
	oldWater += (dt / (l * l)) * fs;
	if (oldWater < 0) {
		oldWater = 0;
	}
	return oldWater;
}

inline void UpdateWaterLevelAndVelocity(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	float oldWater = water[src];
	Flux srcFlux = flux[src];
	Flux neighFlux[4];
	FOR_EACH_DIR_COND(neighs[DIR] != 0, neighFlux[DIR] = flux[neighs[DIR]];)
	float newWater = UpdateWaterLevel(srcFlux, neighs, oldWater, neighFlux);
	float water_level = (oldWater + newWater) * 0.5;
	if (water_level < 0.001) {
		water[src] = newWater;
		velocity[src] = Velocity(0, 0);
		return;
	}
	float dWx = 0, dWy = 0;
	COND_GRID(neighs[0] != 0, dWx -= neighFlux[0].f[2] - srcFlux.f[0];)
	COND_GRID(neighs[1] != 0, dWy -= neighFlux[1].f[3] - srcFlux.f[1];)
	COND_GRID(neighs[2] != 0, dWx += neighFlux[2].f[0] - srcFlux.f[2];)
	COND_GRID(neighs[3] != 0, dWy += neighFlux[3].f[1] - srcFlux.f[3];)
	velocity[src] = Velocity(dWx / (water_level * l), dWy / (water_level * l));
	water[src] = newWater;
}

inline float SinusLocalTiltAngle(int t, int x, int y)
{
	float xl = 2 * l, yl = 2 * l;
	NEIGHBOURS(neighs, x, y);
	if (neighs[0] == 0) {
		neighs[0] = t;
		xl = l;
	} else if (neighs[2] == 0) {
		neighs[2] = t;
		xl = l;
	}
	if (neighs[1] == 0) {
		neighs[1] = t;
		yl = l;
	} else if (neighs[3] == 0) {
		neighs[3] = t;
		yl = l;
	}

	const float dhdx = (TotalGround(neighs[0]) - TotalGround(neighs[2])) / xl;
	const float dhdy = (TotalGround(neighs[1]) - TotalGround(neighs[3])) / yl;
	const float s = dhdx * dhdx + dhdy * dhdy;
	return sqrt(s) / sqrt(1 + s);
}

inline float CalcSedimentCapacity(int src, int x, int y)
{
	Velocity vel = velocity[src];
	const float sinusLocalTiltAngle = SinusLocalTiltAngle(src, x, y);
	const float v = sqrt(vel.x * vel.x + vel.y * vel.y);
	float capacity = Kc * (sinusLocalTiltAngle + v); // * water[src];
	return clamp(capacity, minimumSedimentCapacity, float(1.0));
}

inline void ErosionAndDepositionCalculation(int x, int y)
{
	int src = At(x, y);
	float capacity = CalcSedimentCapacity(src, x, y);
	float sed = sediment[src];

	const float delta = (capacity - sed) * dt;
	float res = 0;
	GroundLayers g = ground[src];

	if (capacity > sed) {
		// picking up sediment
		float f = hardness[1] * delta;
		float l1 = g.layers[1];
		if (f > g.layers[1]) {
			float f2 = (f - l1) * hardness[0] / hardness[1];
			f = g.layers[1] + f2;
		}
		res = f;
	} else {
		// depositing sediment
		res = Kd * delta;
	}
	temp1[src] = res;
}

inline void ErosionAndDepositionUpdate(int x, int y)
{
	int src = At(x, y);
	GroundLayers g = ground[src];
	float sed = sediment[src];
	float ds = temp1[src];
	ground[src] = AddGeneralGround(g, -ds);
	sediment[src] = sed + ds;
}

inline void SedimentTransportation(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	float tmp = 0;
	FOR_EACH_DIR_COND(neighs[DIR] != 0, {
		const float sum = SumFlux(neighs[DIR]);
		if (sum > 0) {
			const float neighSed = sediment[neighs[DIR]];
			const float incomingFlux = flux[neighs[DIR]].f[R_DIR];
			const float dV = dt * neighSed * incomingFlux / sum;
			tmp += dV;
		}
	})
	float sed = sediment[src];
	tmp = tmp + sed * (1.0 - dt);
	temp1[src] = tmp;
}

inline void SedimentTransportationUpdate(int x, int y)
{
	int src = At(x, y);
	sediment[src] = temp1[src];
}

inline void ThermalErosionCalculation(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	NEIGHBOURS_CORNERS(neighCorners, x, y);
	GroundLayers srcGround = ground[src];
	GroundLayers g = srcGround;
	float _h11 = g.layers[1];
	float _h0 = g.layers[0];
	float _h1 = _h0 + _h11;
	float t1 = tangentOfAngleOfRecluse[1] * l;
	float t0 = tangentOfAngleOfRecluse[0] * l;
	float delta = 0;
	const float halfSqrt2 = 0.70710678; // 0.7071067811865475244;

	float n0;
	float n11;
	float n1;

	for (int j = 0; j < 2; ++j) {
		for (int i = 0; i < 4; ++i) {
			if (neighs[i] == 0) {
				continue;
			}

			g = ground[neighs[i]];
			n0 = g.layers[0];
			n11 = g.layers[1];
			n1 = n0 + n11;

			float h0 = _h0;
			float h1 = _h1;
			float h11 = _h11;

			bool swapped = h1 < n1;
			if (swapped) {
				SWAP(n0, h0);
				SWAP(n1, h1);
				SWAP(n11, h11);
			}

			float hh0 = n0 + t0;
			float hh1 = n1 + t1;
			float d = 0;
			if (hh0 < hh1) {
				hh0 = hh1;
			}

			if (h1 > hh1) {
				d = min(h1 - hh1, h11);
				if (h1 > hh0) {
					d = max(d, h1 - hh0);
				}
				delta = delta + (swapped ? d : -d);
			}
		}
		if (j == 1) {
			break;
		}
		for (int i = 0; i < 4; ++i) {
			neighs[i] = neighCorners[i];
		}
		t1 = t1 * 2.0 * halfSqrt2;
		t0 = t0 * 2.0 * halfSqrt2;
	}

	float tmp = delta * (1.0 / 8.0) * 0.5 * dt;
	g = AddGeneralGround(srcGround, tmp);
	velocity[src] = Velocity(g.layers[0], g.layers[1]);
}

inline void ThermalErosionUpdate(int x, int y)
{
	int src = At(x, y);
	Velocity v = velocity[src];
	GroundLayers g;
	g.layers[0] = v.x;
	g.layers[1] = v.y;
	ground[src] = g;
}

inline float EvaporationRate(int x, int y)
{
	return 0.03 * 0.5; // Make it dependent on temperature in place (x,y)
}

inline void Evaporation(int x, int y)
{
		int src = At(x, y);
		temp1[src] = water[src] * (1 - EvaporationRate(x, y) * dt);
}

inline void EvaporationUpdate(int x, int y)
{
		int src = At(x, y);
		water[src] = temp1[src];
}

inline void Smooth(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	const float mult = 128;
	float sum = TotalGround(src) * (mult - 4);
	FOR_EACH_DIR(sum += ((neighs[DIR] != 0) ? TotalGround(neighs[DIR])
											: TotalGround(src));)
	temp1[src] = sum / mult;
}

inline void SmoothUpdate(int x, int y)
{
	int src = At(x, y);
	float dh = temp1[src] - TotalGround(src);
	ground[src].layers[0] += dh;
}

inline void ClearDelta(int x, int y)
{
	int src = At(x, y);
	temp1[src] = 0;
}

#if GL_core_profile
#else
} // namespace HydroPure
#endif
