#if GL_core_profile
#define inline
#define IGNORE_UNUSED(VAR)
#define BUFFER(BLOCK, TYPE, VAR)                                               \
	buffer BLOCK { TYPE VAR[]; };
#define static
#else
#include <cmath>
#define layout(a)
#define uniform
#define IGNORE_UNUSED(VAR) (void)VAR;
#define BUFFER(BLOCK, TYPE, VAR) static TYPE *VAR = nullptr;
inline float min(float A, float B) { return A < B ? A : B; }
inline float max(float A, float B) { return A > B ? A : B; }
#endif

#define COND_GRID(COND, EXPR)                                                  \
	if (COND) {                                                                \
		EXPR                                                                   \
	}

#define NEIGHBOURS(VAR_NAME, X, Y)                                             \
	int VAR_NAME[4] = {Neighbour(X, Y, 0), Neighbour(X, Y, 1),                 \
					   Neighbour(X, Y, 2), Neighbour(X, Y, 3)}

#define NEIGHBOURS_CORNERS(VAR_NAME, X, Y)                                     \
	int VAR_NAME[4] = {                                                        \
		At(X - 1, Y - 1),                                                      \
		At(X + 1, Y - 1),                                                      \
		At(X - 1, Y + 1),                                                      \
		At(X + 1, Y + 1),                                                      \
	}

#define FOR_EACH_DIR(CODE)                                                     \
	{                                                                          \
		for (int DIR = 0; DIR < 4; ++DIR) {                                    \
			int R_DIR = (DIR + 2) & 3;                                         \
			IGNORE_UNUSED(R_DIR)                                               \
			CODE                                                               \
		}                                                                      \
	}

#define FOR_EACH_DIR_COND(COND, CODE) FOR_EACH_DIR(COND_GRID(COND, CODE))

#define REV_DIR(DIR) ((DIR + 2) & 3)

#define SWAP(A, B)                                                             \
	{                                                                          \
		float tmp = A;                                                         \
		A = B;                                                                 \
		B = tmp;                                                               \
	}

#if GL_core_profile
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
#else
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
namespace HydroPure
{
#endif

layout(binding = 1) BUFFER(BufferBlock1, GroundLayers, ground);
layout(binding = 2) BUFFER(BufferBlock2, float, water);
layout(binding = 3) BUFFER(BufferBlock3, float, sediment);
layout(binding = 4) BUFFER(BufferBlock4, float, temp1);
layout(binding = 5) BUFFER(BufferBlock5, Velocity, velocity);
layout(binding = 6) BUFFER(BufferBlock6, Flux, flux);

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

inline void AddGeneralGround(int t, float dv)
{
	GroundLayers g = ground[t];
	g.layers[1] += dv;
	if (g.layers[1] < 0) {
		g.layers[0] += g.layers[1];
		g.layers[1] = 0;
	}
	ground[t] = g;
}

inline int At(int x, int y)
{
	if (x < 0)
		return 0; // At(x+width, y); // 0;
	else if (x >= width)
		return 0; // At(x-width, y); // 0;
	if (y < 0)
		return 0; // At(x, y+height); // 0;
	else if (y >= height)
		return 0; // At(x, y-height); // 0;
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
	// TODO: should not happen
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

inline Flux LimitFlux(int src, Flux f)
{
	const float outflux = SumFluxF(f);
	const float _water = water[src] * l * l;
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
	NEIGHBOURS(neighs, x, y);

	float srcSum = TotalGround(src) + sediment[src] + water[src];

	FOR_EACH_DIR_COND(
		neighs[DIR] != 0,
		(f.f[DIR] = CalcFluxInDirection(src, neighs[DIR], DIR, f, srcSum));)
	flux[src] = LimitFlux(src, f);
}

inline float UpdateWaterLevel(int src, int neighs[4], float oldWater)
{
	float fs = 0;
	FOR_EACH_DIR_COND(neighs[DIR] != 0,
					  fs += flux[neighs[DIR]].f[R_DIR] - flux[src].f[DIR];)
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
	float water_level = water[src];
	float newWater = UpdateWaterLevel(src, neighs, water_level);
	water_level = (water_level + newWater) * 0.5;
	if (water_level < 0.001) {
		velocity[src].x = 0;
		velocity[src].y = 0;
		return;
	}
	float dWx = 0, dWy = 0;
	COND_GRID(neighs[0] != 0, dWx -= flux[neighs[0]].f[2] - flux[src].f[0];)
	COND_GRID(neighs[1] != 0, dWy -= flux[neighs[1]].f[3] - flux[src].f[1];)
	COND_GRID(neighs[2] != 0, dWx += flux[neighs[2]].f[0] - flux[src].f[2];)
	COND_GRID(neighs[3] != 0, dWy += flux[neighs[3]].f[1] - flux[src].f[3];)
	velocity[src].x = dWx / (water_level * l);
	velocity[src].y = dWy / (water_level * l);
	water[src] = newWater;
}

inline float SinusLocalTiltAngle(int t, int x, int y)
{
	float xl = 2 * l, yl = 2 * l;
	NEIGHBOURS(neighs, x, y);
	if (neighs[0] == 0) {
		neighs[0] = t;
		xl = l;
	}
	if (neighs[1] == 0) {
		neighs[1] = t;
		yl = l;
	}
	if (neighs[2] == 0) {
		neighs[2] = t;
		xl = l;
	}
	if (neighs[3] == 0) {
		neighs[3] = t;
		yl = l;
	}

	const float dhdx = (TotalGround(neighs[0]) - TotalGround(neighs[2])) / xl;
	const float dhdy = (TotalGround(neighs[1]) - TotalGround(neighs[3])) / yl;
	const float s = dhdx * dhdx + dhdy * dhdy;
	return sqrt(s) / sqrt(1 + s);
}

inline void ErosionAndDepositionCalculation(int x, int y)
{
	int src = At(x, y);
	Velocity vel = velocity[src];
	const float sinusLocalTiltAngle = SinusLocalTiltAngle(src, x, y);
	const float v = sqrt(vel.x * vel.x + vel.y * vel.y);
	float capacity = Kc * sinusLocalTiltAngle * v;
	if (capacity < minimumSedimentCapacity)
		capacity = minimumSedimentCapacity;
	if (capacity > 1)
		capacity = 1;

	const float delta = (capacity - sediment[src]) * dt;
	float res = 0;
	GroundLayers g = ground[src];

	if (capacity > sediment[src]) {
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
	float ds = temp1[src];
	AddGeneralGround(src, -ds);
	sediment[src] += ds;
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
	// 	float f = SumFlux(src);
	// 	float sed = sediment[src];
	// 	if (f > 0) {
	// 		sed *= 1 - dt;
	// 	}
	sediment[src] = /*sed */ +temp1[src];
}

inline void ThermalErosionCalculation(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	NEIGHBOURS_CORNERS(neighCorners, x, y);
	GroundLayers g = ground[src];
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
	temp1[src] = delta * (1.0 / 8.0) * 0.5 * dt;
}

inline void ThermalErosionUpdate(int x, int y)
{
	int src = At(x, y);
	AddGeneralGround(src, temp1[src]);
}

inline float EvaporationRate(int x, int y)
{
	return 0.04; // Make it dependent on temperature in place (x,y)
}

inline void Evaporation(int x, int y)
{
	int src = At(x, y);
	water[src] *= (1 - EvaporationRate(x, y) * dt);
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
