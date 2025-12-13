
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
%:include <atomic>
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
layout(std430, binding = 4) buffer BufferBlock71
{
	layout(offset = (PADDING * 1 + ELEMENTS * 0) * 4) float water[];
};
layout(std430, binding = 4) buffer BufferBlock72
{
	layout(offset = (PADDING * 2 + ELEMENTS * 1) * 4) float sediment[];
};
layout(std430, binding = 4) buffer BufferBlock72_2
{
	layout(offset = (PADDING * 2 + ELEMENTS * 1) * 4) int sediment_int[];
};
layout(std430, binding = 4) buffer BufferBlock73
{
	layout(offset = (PADDING * 3 + ELEMENTS * 2) * 4) float temp1[];
};
layout(std430, binding = 4) buffer BufferBlock73_2
{
	layout(offset = (PADDING * 3 + ELEMENTS * 2) * 4) int temp1_int[];
};
layout(std430, binding = 4) buffer BufferBlock74
{
	layout(offset = (PADDING * 4 + ELEMENTS * 3) * 4) float temp2[];
};
layout(std430, binding = 4) buffer BufferBlock74_2
{
	layout(offset = (PADDING * 4 + ELEMENTS * 3) * 4) int temp2_int[];
};
#else
static float *water = nullptr;
static float *sediment = nullptr;
static int *sediment_int = nullptr;
static float *temp1 = nullptr;
static int *temp1_int = nullptr;
static float *temp2 = nullptr;
static int *temp2_int = nullptr;
#endif

#if GL_core_profile
#else
namespace _____holder
{
inline void _holder()
{
	(void)iteration;
	(void)sediment_int;
	(void)temp1_int;
	(void)temp2;
	(void)temp2_int;
}
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

inline int AtClamp(int x, int y)
{
	x = clamp(x, 0, width-1);
	y = clamp(y, 0, height-1);
	return (x * height + y) + 1;
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

#if GL_core_profile
#else
inline int atomicCompSwap(int &mem, int cmp, int val)
{
	static_assert(sizeof(int) == sizeof(std::atomic<int>));
	std::atomic<int> *atom = (std::atomic<int> *)&mem;
	atom->compare_exchange_weak(cmp, val);
	return cmp;
}

inline int atomicExchange(int &mem, int val)
{
	static_assert(sizeof(int) == sizeof(std::atomic<int>));
	std::atomic<int> *atom = (std::atomic<int> *)&mem;
	return atom->exchange(val);
}

inline float intBitsToFloat(int val)
{
	return *(float *)&val;
}

inline int floatBitsToInt(float val)
{
	return *(int *)&val;
}
#endif

#define GetAtomicFloat(MEM) intBitsToFloat(atomicExchange(MEM, 0))
#define SetAtomicFloat(MEM, VAL) atomicExchange(MEM, floatBitsToInt(VAL))

// returns 1 when lock failed
// returns 0 when lock acquired
inline int TryLock(int id)
{
	int old = atomicCompSwap(temp2_int[id], 0, 1);
	return old;
}

inline void Unlock(int id) { atomicExchange(temp2_int[id], 0); }

#define LOCK_LOOP(ID, CODE)                                                    \
	{                                                                          \
		int lockAvailable = 0;                                                 \
		do {                                                                   \
			lockAvailable = TryLock(ID);                                       \
			if (lockAvailable == 0) {                                          \
				CODE Unlock(ID);                                               \
			}                                                                  \
		} while (lockAvailable != 0);                                          \
	}

struct RiverSource {
	ivec2 coord;
	float amount;
};

inline RiverSource CalcRiverSource(ivec2 p, uint gridSize, uint seed)
{
	const uint sourcesGridSize = gridSize;
	uvec2 v = uvec2(p);
	uvec2 md = v % sourcesGridSize;
	uvec2 base = v / sourcesGridSize;

	uvec4 rnd = RandomUint(ivec3(int(base.x), int(base.y), seed));

	uvec2 source = v - md + (uvec2(rnd.x, rnd.y) % sourcesGridSize);
	if (rnd.z % 7 == 0) {
		return RiverSource(ivec2(source), 0.0);
	}
	float amount = float(rnd.w % 1511u) / 312.0;
// 	amount = amount * amount * amount;
// 	amount *= amount;
	amount += 1.0;
	return RiverSource(ivec2(source), amount);
}

inline float CalcRain(int x, int y, int src)
{
	ivec2 coord = ivec2(x, y);
	RiverSource s = CalcRiverSource(coord, 17, iteration);
	if (length(vec2(s.coord-coord)) < 5) {
		return 0.0001;
	}
	return 0;
	
	
	
	
	float val;

	// 	uvec4 rnd = RandomUint(ivec3(x, y, iteration));
	// 	val = (rnd.x % 321381) / 321380.0;

	val = sin(iteration / 5.0 * dt);
	val = clamp(val, float(0.0), float(1.0));

	/*
	vec3 coord =
		vec3(float(x), float(y), float(iteration % 100000) * 0.02) *
	float(0.01); vec3 grad; val = (psrdnoise(coord, vec3(0, 0, 0), 0.001 *
	float(iteration % 100000), grad) * 0.5 + 0.5); val = clamp(val, float(0.0),
	float(1.0));
	*/

	val *= val;
	val *= val;
	val *= val;
	val *= 0.002;
	val *= 4.0;
	return val;
}

inline void RainAndRiverUpdate(int x, int y)
{
	int src = At(x, y);
	float w = water[src];

	float rain = CalcRain(x, y, src);
	float sources = 0.0;

	RiverSource rs = CalcRiverSource(ivec2(x, y), 73, 43124324);
	ivec2 ap = (rs.coord - ivec2(x, y));
	int dp = ap.x * ap.x + ap.y * ap.y;
	if (dp <= 1) {
		sources = rs.amount / (1.0 + dp); // * (1.0 / 5.0);
	}

	water[src] = w + (rain + sources) * dt;
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
	f = LimitFlux(src, f, waterLevel);
	flux[src] = f;
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

inline void WaterLevelAndVelocityUpdate(int x, int y)
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

inline vec3 Normal(int x, int y)
{
	const float halfSqrt2 = 0.70710678; // 0.7071067811865475244;
	const float sqrt2 = 2.0 * halfSqrt2;
	float h[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	float g;
	for (int i=0; i<3; ++i) {
		for (int j=0; j<3; ++j) {
			int n = AtClamp(i+x-1, j+y-1);
			h[i][j] = TotalGround(n);
		}
	}

	vec3 n1 = vec3(0, 0, 0);
	g = h[1][1] - h[2][1]; n1 += normalize(vec3(g, l, 0));
	g = h[0][1] - h[1][1]; n1 += normalize(vec3(g, l, 0));
	g = h[1][1] - h[1][2]; n1 += normalize(vec3(0, l, g));
	g = h[1][0] - h[1][1]; n1 += normalize(vec3(0, l, g));
	n1 *= 0.15;

	vec3 n2 = vec3(0, 0, 0);
	g = h[0][0] - h[1][1]; n2 += normalize(vec3(+g, l*sqrt2, +g));
	g = h[2][0] - h[1][1]; n2 += normalize(vec3(-g, l*sqrt2, +g));
	g = h[0][2] - h[1][1]; n2 += normalize(vec3(+g, l*sqrt2, -g));
	g = h[2][2] - h[1][1]; n2 += normalize(vec3(-g, l*sqrt2, -g));
	n2 *= 0.1;

	return n1 + n2;
}

inline vec3 NormalG(int x, int y)
{
	int t = At(x, y);
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

// 	const float dhdx = (TotalGround(neighs[0]) - TotalGround(neighs[2])) / xl;
// 	const float dhdy = (TotalGround(neighs[1]) - TotalGround(neighs[3])) / yl;
// 	const float s = dhdx * dhdx + dhdy * dhdy;
// 	return sqrt(s) / sqrt(1 + s);

// 	return (
// 	normalize(vec3((TotalGround(neighs[0]) - TotalGround(neighs[2])), xl, 0)) +
// 	normalize(vec3(0, yl, (TotalGround(neighs[1]) - TotalGround(neighs[3]))))
// 	) * float(0.5);
	
	const vec3 ax = vec3(xl, (TotalGround(neighs[0]) - TotalGround(neighs[2])), 0);
	const vec3 ay = vec3(0, (TotalGround(neighs[1]) - TotalGround(neighs[3])), yl);
	vec3 n = cross(ax, ay);
	if (n.y < 0) {
		n = -n;
	}
	return n;
}

inline float SinusLocalTiltAngle(int t, int x, int y)
{
	vec3 cr = Normal(x, y);
	const float _dot = cr.y; // dot(cr, vec3(0,1,0));
	const float _cos = _dot / length(cr);
	return _cos;
}

inline float CalcSedimentCapacity(int src, int x, int y)
{
	Velocity vel = velocity[src];
	float sinusLocalTiltAngle =
// 		clamp(
				SinusLocalTiltAngle(src, x, y)
// 				, float(0), float(1))
		;
	if (
			sinusLocalTiltAngle >= 0
// 			&&
// 			sinusLocalTiltAngle <= 1.01
			) {
	} else {
		ground[src].layers[0] = -1000;
	}
	float w = clamp(water[src], float(0), float(1));
	const float v =
	//	clamp(
			sqrt(vel.x * vel.x + vel.y * vel.y)
	//		, minimumSedimentCapacity, float(10))
		;
	float capacity = Kc * sinusLocalTiltAngle * v * w;
	return capacity;
// 	return clamp(capacity, minimumSedimentCapacity, float(100.0));
	return clamp(capacity, minimumSedimentCapacity, float(1.0));
// 	return clamp(capacity, minimumSedimentCapacity, w * Kc);
}

inline void ErosionAndDepositionCalculation(int x, int y)
{
	int src = At(x, y);
	float capacity = CalcSedimentCapacity(src, x, y);
	float sed = sediment[src];
	temp2_int[src] = 0;

	const float delta = (capacity - sed) * dt;
	float res = 0;
	GroundLayers g = ground[src];

	if (capacity > sed) {
		// picking up sediment
		float f = hardness[1] * delta / dt;
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
	g = AddGeneralGround(g, -res);
	sediment[src] = sed + res;
	temp1[src] = g.layers[0];
	temp2[src] = g.layers[1];
}

inline void ErosionAndDepositionUpdate(int x, int y)
{
	int src = At(x, y);
	GroundLayers g;
	g.layers[0] = temp1[src];
	g.layers[1] = temp2[src];
	ground[src] = g;
	temp1[src] = 0;
	temp2_int[src] = 0;
}

#define R 2
#define D (R * 2 + 1)

inline vec2 VelocityDx(Velocity vel)
{
	const Velocity _v = vel;
	const vec2 v1 = -vec2(_v.x, _v.y) * dt;
	const float len = length(v1);
	if (len > D) {
		return v1 * (float(D) / len);
	}
	return v1;
}

inline void SedimentTransportation(int x, int y)
{
	int src = At(x, y);


	/*
	NEIGHBOURS(neighs, x, y);
	float tmp = 0;
	FOR_EACH_DIR_COND(neighs[DIR] != 0, {
		const float sum = SumFlux(neighs[DIR]);
		if (sum > 0) {
			const float neighSed = sediment[neighs[DIR]];
			const float incomingFlux = flux[neighs[DIR]].f[R_DIR];
			const float dV = neighSed * incomingFlux / sum;
			tmp += dV;
		}
	})
	const float sed = sediment[src];
// 	tmp = tmp * dt + sed * (1.0 - dt);
	tmp = sed + (tmp - sed) * dt;
	temp1[src] = tmp;
	return;
	*/

	
	/*
	float _sed[D][D];
	vec2 srcVel;
	float deltaSed = 0.0;
	for (int i = 0; i < D; ++i) {
		for (int j = 0; j < D; ++j) {
			int id = At(i + x - R, j + y - R);
			if (id != 0) {
				_sed[i][j] = sediment[id];
				vec2 v = VelocityDx(velocity[id]);
				if (i == j && j == R) {
					srcVel = v;
				}
				vec2 _p = vec2(i, j) + v;
				vec2 dist = abs(_p - vec2(R, R));
				if (dist.x < 1.0 && dist.y < 1.0) {
					dist = vec2(1, 1) - dist;
					float fr = dist.x * dist.y;
					deltaSed -= fr;
				}
			} else {
				_sed[i][j] = 0;
			}
		}
	}
	deltaSed *= _sed[R][R];

	{
		const vec2 p = srcVel + vec2(R, R);
		const ivec2 ip = ivec2(floor(p));
		const vec2 f = p - vec2(ip);
		
		deltaSed += _sed[ip.x][ip.y] * (1 - f.x) * (1 - f.y);
		deltaSed += _sed[ip.x + 1][ip.y] * (f.x) * (1 - f.y);
		deltaSed += _sed[ip.x][ip.y + 1] * (1 - f.x) * (f.y);
		deltaSed += _sed[ip.x + 1][ip.y + 1] * (f.x) * (f.y);
	}

	const float tmp = _sed[R][R] + deltaSed / 9.0;
	temp1[src] = tmp;
	*/

	
	
	/*
	const Velocity v = velocity[src];
	const vec2 srcVel = -vec2(v.x, v.y) * dt;
	float deltaSed = 0.0;
	
	const vec2 p = srcVel + vec2(x, y);
	const ivec2 ip = ivec2(floor(p));
	const vec2 f = p - vec2(ip);
	
	const ivec2 is[4] = {{0,0}, {1,0}, {0,1}, {1,1}};
	for (int i=0; i<4; ++i) {
		const ivec2 pi = is[i] + ip;//ivec2(x, y);
		const int id = At(pi.x, pi.y);
		if (id != 0) {
			const float fx = is[i].x == 0 ? float(1) - f.x : f.x;
			const float fy = is[i].y == 0 ? float(1) - f.y : f.y;
			deltaSed += sediment[id] * fx * fy;
		}
	}
	temp1[src] = (sediment[src] + deltaSed) * 0.5;
	*/
	

	/*
	const Velocity v = velocity[src];
	const vec2 srcVel = -vec2(v.x, v.y) * dt;
	float deltaSed = 0.0;
	
	const vec2 p = srcVel + vec2(x, y);
	const ivec2 ip = ivec2(floor(p));
	const vec2 f = p - vec2(ip);

	const int ids[4] = {
		At(ip.x, ip.y),
		At(ip.x+1, ip.y),
		At(ip.x, ip.y+1),
		At(ip.x+1, ip.y+1)};

	if (ids[0] != 0) deltaSed += sediment[ids[0]] * (1 - f.x) * (1 - f.y);
	if (ids[1] != 0) deltaSed += sediment[ids[1]] * (f.x) * (1 - f.y);
	if (ids[2] != 0) deltaSed += sediment[ids[2]] * (1 - f.x) * (f.y);
	if (ids[3] != 0) deltaSed += sediment[ids[3]] * (f.x) * (f.y);
	temp1[src] = deltaSed;
	*/


	/*
	const Velocity v = velocity[src];
	const vec2 srcVel = -vec2(v.x, v.y) * dt;
	float deltaSed = 0.0;
	
	const vec2 p = srcVel + vec2(x, y);
	const ivec2 ip = ivec2(floor(p));
	const vec2 f = p - vec2(ip);

	const int ids[4] = {
		At(ip.x, ip.y),
		At(ip.x+1, ip.y),
		At(ip.x, ip.y+1),
		At(ip.x+1, ip.y+1)};
	
	const float factors[4] = {
		(1 - f.x) * (1 - f.y),
		(f.x) * (1 - f.y),
		(1 - f.x) * (f.y),
		(f.x) * (f.y)
	};

	for (int i = 0; i < 4; ++i) {
		int id = ids[i];
		float f = factors[i];
		if (f >= 0 && f <= 1) {
		} else {
			ground[src].layers[0] = -10000.0;
		}
		if (id != 0) {
			float ds = 0;
			LOCK_LOOP(id,
				{
					float s = GetAtomicFloat(sediment_int[id]);
					if (s < 0) {
						ds = 0;
						SetAtomicFloat(sediment_int[id], 0.0);
					} else {
						ds = s * f;
						SetAtomicFloat(sediment_int[id], s - ds);
					}
				});
			deltaSed += ds;
		}
	}
	temp1[src] = deltaSed;
	*/


	/*
	if (false)
	{
		const Velocity v = velocity[src];
		const vec2 srcVel = vec2(v.x, v.y) * dt;
		float sed = sediment[src];

		const vec2 p = srcVel + vec2(x, y);
		const ivec2 ip = ivec2(floor(p));
		const vec2 f = p - vec2(ip);

		const int ids[4] = {
			At(ip.x, ip.y),
			At(ip.x+1, ip.y),
			At(ip.x, ip.y+1),
			At(ip.x+1, ip.y+1)};

		const float factors[4] = {
			(1 - f.x) * (1 - f.y),
			(f.x) * (1 - f.y),
			(1 - f.x) * (f.y),
			(f.x) * (f.y)
		};

		for (int i = 0; i < 4; ++i) {
			int id = ids[i];
			float f = factors[i];
			if (f >= 0 && f <= 1) {
			} else {
				ground[src].layers[0] = -10000.0;
			}
			if (id != 0) {
				float ds = 0;
				LOCK_LOOP(id,
					{
						float s = GetAtomicFloat(temp1_int[id]);
						s += f * sed;
						if (s < 0) {
							ds = 0;
							SetAtomicFloat(temp1_int[id], 0.0);
						} else {
							ds = s * f;
							SetAtomicFloat(temp1_int[id], s - ds);
						}
					});
			}
		}
		return;
	}
	*/



	const Velocity v = velocity[src];
	const vec2 srcVel = -vec2(v.x, v.y) * dt;
	const vec2 p = srcVel + vec2(x, y);
	const ivec2 ip = ivec2(floor(p));
	const vec2 f = p - vec2(ip);

	const int ids[4] = {
		At(ip.x, ip.y),
		At(ip.x+1, ip.y),
		At(ip.x, ip.y+1),
		At(ip.x+1, ip.y+1)};
	const float factors[4] = {
		(1 - f.x) * (1 - f.y),
		(f.x) * (1 - f.y),
		(1 - f.x) * (f.y),
		(f.x) * (f.y)
	};

	for (int i = 0; i < 4; ++i) {
		int id = ids[i];
		float f = factors[i];
		if (id != 0) {
			float ds = 0;
			float othSed = sediment[id];
			ds = othSed * f;
			LOCK_LOOP(id,
				{
// 					float s = temp1[id];
// 					s += ds;
// 					temp1[id] = s;
					float s = GetAtomicFloat(temp1_int[id]);
					s += ds;
					SetAtomicFloat(temp1_int[id], s);
				});
		}
	}
}
#undef R
#undef D

inline void SedimentTransportation_Stage2(int x, int y)
{
	int src = At(x, y);
	const Velocity v = velocity[src];
	const vec2 srcVel = -vec2(v.x, v.y) * dt;
	float deltaSed = 0.0;
	const vec2 p = srcVel + vec2(x, y);
	const ivec2 ip = ivec2(floor(p));
	const vec2 f = p - vec2(ip);

	const int ids[4] = {
		At(ip.x, ip.y),
		At(ip.x+1, ip.y),
		At(ip.x, ip.y+1),
		At(ip.x+1, ip.y+1)};
	const float factors[4] = {
		(1 - f.x) * (1 - f.y),
		(f.x) * (1 - f.y),
		(1 - f.x) * (f.y),
		(f.x) * (f.y)
	};

	for (int i = 0; i < 4; ++i) {
		int id = ids[i];
		float f = factors[i];
		if (id != 0) {
			float ds = 0;
			float othSed = sediment[id];
			ds = othSed * f;
			float othSum = temp1[id];
			
			if (othSum > 0.0) {
				deltaSed += ds / (othSum > othSed ? othSum : 1.0);
			}
		}
	}

	float sedSrc = sediment[src];
	float sedSrcSum = temp1[src];
	if (sedSrcSum > 0.0) {
		temp2[src] = deltaSed;
	} else {
		temp2[src] = sedSrc + deltaSed;
	}
}

inline void SedimentTransportationUpdate(int x, int y)
{
	int src = At(x, y);
	sediment[src] = temp2[src];
// 	sediment[src] = temp1[src];
// 	sediment[src] += temp1[src];
}

/*
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

	float tmp = delta * (1.0 / (4.0 + 8.0 * halfSqrt2)) * 0.5;// * dt;
	g = AddGeneralGround(srcGround, tmp);
	velocity[src] = Velocity(g.layers[0], g.layers[1]);
}
*/

inline ivec2 mul(const ivec2 mat[2], ivec2 vec)
{
	return ivec2(dot(vec, mat[0]), dot(vec, mat[1]));
}

inline int GetThermalRadius()
{
// 	if (iteration > 1000) {
// 		return 0;
// 	}
// 	return 5;
	
	return ((iteration%10) == 0 || iteration < 10000) ? 5 : 0;
	
	if (iteration < 10000) {
		return 5;
	} else if (iteration < 1100 || (iteration % 2301 < 6)) {
		return 10;
	}
	return 0;
	
	int TER = 4;
	if (iteration < 100 || (iteration % 2301 < 6)) {
		TER = 50;
	} else if (iteration < 100) {
		TER = 25;
	} else if (iteration < 200) {
		TER = 20;
	} else if (iteration % 1000 == 7) {
		TER = 15;
	} else if (iteration % 133 == 23) {
		TER = 8;
	} else if (iteration % 17 == 9) {
		TER = 6;
	} else if (iteration % 7 == 4) {
		TER = 5;
	}
	if (TER < 50) {
		return 0;
	}
	return TER;
}

inline float ThermalErosionWaterFactorFunction(float a, float b)
{
	const float x = (a + b) / 10;
	return 1.0 / (1.0 + x);
// 	return 1.0 / (1.0 + x * x);
// 	return 0.5 / (1.0 + x * x) + 0.5;
}

inline void ThermalErosionCalculation(int x, int y)
{
	int src = At(x, y);

	GroundLayers srcGround = ground[src];
	GroundLayers g = srcGround;
	const float srcWater = water[src];

	int TER = GetThermalRadius();
	if (TER <= 0) {
		return;
	}

	float _h11 = g.layers[1];
	float _h0 = g.layers[0];
	float _h1 = _h0 + _h11;
	float _t1 = tangentOfAngleOfRecluse[1] * l;
	float _t0 = tangentOfAngleOfRecluse[0] * l;
	float delta = 0;

	const ivec2 rots[4][2] = {{{1, 0}, {0, 1}},
							  {{0, -1}, {1, 0}},
							  {{-1, 0}, {0, -1}},
							  {{0, 1}, {-1, 0}}};
	float sum = 0.0;

	for (int _i = 0; _i <= TER; ++_i) {
		for (int _j = 1; _j <= TER; ++_j) {
			ivec2 _pi = ivec2(_i, _j);
			float dist = sqrt(_i * _i + _j * _j);
			if (dist + 1 >= TER) {
				continue;
			}
			sum += 4;
// 			sum += dist * 4.0 * l;
			float t0 = dist * _t0;
			float t1 = dist * _t1;
			for (int k = 0; k < 4; ++k) {
				ivec2 c = mul(rots[k], _pi);
				ivec2 p = c + ivec2(x, y);

				int neigh = At(p.x, p.y);
				if (neigh == 0 || neigh == src) {
					continue;
				}
				
				const float nw = water[neigh];

				g = ground[neigh];
				float n0 = g.layers[0];
				float n11 = g.layers[1];
				float n1 = n0 + n11;

				float h0 = _h0;
				float h1 = _h1;
				float h11 = _h11;

				bool swapped = h1 < n1;
				if (swapped) {
					SWAP(n0, h0);
					SWAP(n11, h11);
					SWAP(n1, h1);
				}

				float f = ThermalErosionWaterFactorFunction(srcWater, nw);
				float hh0 = n0 + t0;
				float hh1 = n1 + t1 * f;
				float d = 0;
				if (hh0 < hh1) {
					hh0 = hh1;
				}

				if (h1 > hh1) {
					d = min(h1 - hh1, h11);
					if (h1 > hh0) {
						d = max(d, h1 - hh0);
					}
					delta += (swapped ? d : -d); // (TER < 20 ? dist : 1.0);
				}
			}
		}
	}

// 	sum = (TER * 2 + 1) * (TER * 2 + 1) - 1;
	float tmp = delta * (1.0 / sum) * 0.5;// * (iteration > 10000 ? 0.2 : 1.0); // * dt;
	g = AddGeneralGround(srcGround, tmp);
	velocity[src] = Velocity(g.layers[0], g.layers[1]);
}

inline void ThermalErosionUpdate(int x, int y)
{
	int TER = GetThermalRadius();
	if (TER <= 0) {
		return;
	}

	int src = At(x, y);
	Velocity v = velocity[src];
	GroundLayers g;
	g.layers[0] = v.x;
	g.layers[1] = v.y;
	ground[src] = g;
}

inline float EvaporationRate(int x, int y)
{
	return 0.03 / 16.0; // / 16.0; // Make it dependent on temperature in place (x,y)
}

inline void Evaporation(int x, int y)
{
	int src = At(x, y);
	float w = water[src];
	float evap = EvaporationRate(x, y) * dt;
	evap = clamp(evap, float(0.0), w);
	water[src] = w - evap;
}

// TODO: replace with selectional smoothing, to smooth only where slope
// changes very rapidly
inline void Smooth(int x, int y)
{
	int src = At(x, y);
	NEIGHBOURS(neighs, x, y);
	NEIGHBOURS_CORNERS(neighCorners, x, y);
	const float mult = 256;
	const float halfSqrt2 = 0.70710678; // 0.7071067811865475244;

	float tg = TotalGround(src);
	float sum = tg * (mult - 4 - 4.0 * halfSqrt2);
	float hs = 0;
	FOR_EACH_DIR(hs += ((neighs[DIR] != 0) ? TotalGround(neighs[DIR]) : tg);)
	FOR_EACH_DIR(
		hs += ((neighCorners[DIR] != 0) ? TotalGround(neighCorners[DIR]) : tg) *
			  halfSqrt2;)
	sum += hs;
	temp1[src] = ground[src].layers[0] + sum / mult - tg;
}

inline void SmoothUpdate(int x, int y)
{
	int src = At(x, y);
	float dh = temp1[src]; // - TotalGround(src);
	ground[src].layers[0] = dh;
}

#if GL_core_profile
#else
} // namespace HydroPure
#endif
