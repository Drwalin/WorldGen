// #include <cmath>
#include <cassert>

#include <chrono>
#include <vector>
#include <thread>
#include <atomic>
#include <functional>

#include "../include/worldgen/HydroErosion.hpp"
// #include "../include/worldgen/HydroErosionMacros.hpp"

#include "HydroErosionPure.hpp"

std::atomic<float> f;

int a = (f += 1.31, 0);


/*
inline float Grid::SumFlux(int t) {
	return flux[t].f[0] + flux[t].f[1] + flux[t].f[2] + flux[t].f[3];
}

template<bool safe, int dir>
inline float Grid::CalcFluxInDirection(int src, int neigh) const {
	if(neigh == 0)
		return 0.0f;
	int dst = neigh;
	float dh = (ground[src].Total() - ground[dst].Total()) + (sediment[src] - sediment[dst]) + (water[src] - water[dst]);
	float f = flux[src].fluxArray[dir] + dt * A * g * dh / l;
	if(f < 0.0f)
		f = 0.0f;
	return f;
}

inline void Grid::LimitFlux(int src) {
	const float outflux = SumFlux(src);
	const float water = this->water[src] * l*l;
	if(outflux <= 0.001) {
		FOR_EACH_DIR({flux[src].fluxArray[DIR] = 0.0f;});
		return;
	}
	float K = water / (outflux * dt);
	if (K > 1.0f) {
		K = 1.0f;
	}
	if (K >= 0.0f) {
		FOR_EACH_DIR({flux[src].fluxArray[DIR] *= K;});
	}
}
*/

template<bool safe>
void Grid::CalcOutflux(int x, int y) {
	HydroPure::CalcOutflux(x, y);
	/*
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], (flux[src].fluxArray[DIR] = CalcFluxInDirection<safe, DIR>(src, neighs[DIR])));
	LimitFlux(src);
	*/
}

/*
template<bool safe>
void Grid::UpdateWaterLevel(int src, int* neighs) {
	float fs = 0;
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], fs += flux[neighs[DIR]].fluxArray[R_DIR] - flux[src].fluxArray[DIR]);
	water[src] += (dt/(l*l)) * fs;
	if (water[src] < 0.0f) {
		water[src] = 0;
	}
}
*/

template<bool safe>
void Grid::UpdateWaterLevelAndVelocity(int x, int y) {
	HydroPure::UpdateWaterLevelAndVelocity(x, y);
	/*
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	float water_level = water[src];
	UpdateWaterLevel<safe>(src, neighs);
	water_level = (water_level + water[src]) * 0.5f;
	if(water_level < 0.001) {
		velocity[src].x = 0;
		velocity[src].y = 0;
		return;
	}
	float dWx = 0.0f, dWy = 0.0f;
	SAFE_COND_GRID(neighs[0], dWx -= flux[neighs[0]].fluxArray[2] - flux[src].fluxArray[0]);
	SAFE_COND_GRID(neighs[1], dWy -= flux[neighs[1]].fluxArray[3] - flux[src].fluxArray[1]);
	SAFE_COND_GRID(neighs[2], dWx += flux[neighs[2]].fluxArray[0] - flux[src].fluxArray[2]);
	SAFE_COND_GRID(neighs[3], dWy += flux[neighs[3]].fluxArray[1] - flux[src].fluxArray[3]);
	velocity[src].x = dWx / (water_level * l);
	velocity[src].y = dWy / (water_level * l);
	*/
}

/*
template<bool safe>
inline float Grid::SinusLocalTiltAngle(int t, int x, int y) {
	float xl, yl;
	NEIGHBOURS(neighs, x, y);
	if constexpr (safe) {
		xl = yl = 2.0f * l;
		if(neighs[0] == 0) { neighs[0] = t; xl = l; }
		if(neighs[1] == 0) { neighs[1] = t; yl = l; }
		if(neighs[2] == 0) { neighs[2] = t; xl = l; }
		if(neighs[3] == 0) { neighs[3] = t; yl = l; }
	} else {
		xl = yl = l;
	}
	
	const float dhdx = (ground[neighs[0]].Total() - ground[neighs[2]].Total()) / xl;
	const float dhdy = (ground[neighs[1]].Total() - ground[neighs[3]].Total()) / yl;
	const float s = dhdx*dhdx + dhdy*dhdy;
	return sqrt(s) / sqrt(1.0f + s);
}
*/

template<bool safe>
inline void Grid::ErosionAndDepositionCalculation(int x, int y) {
	HydroPure::ErosionAndDepositionCalculation(x, y);
	/*
	int src = At<false>(x, y);
	const float sinusLocalTiltAngle = SinusLocalTiltAngle<safe>(src, x, y);
	const float v = sqrt(velocity[src].x * velocity[src].x + velocity[src].y * velocity[src].y);
	float capacity = Kc * sinusLocalTiltAngle * v;
	if(capacity < minimumSedimentCapacity)
		capacity = minimumSedimentCapacity;
	if (capacity > 1.0f)
		capacity = 1.0f;
	
	const float delta = (capacity - sediment[src]) * dt;
	
	if(capacity > sediment[src]) {
		// picking up sediment
		float f = hardness[1] * delta;
		float l1 = ground[src][1];
		if (f > ground[src][1]) {
			float f2 = (f - l1) * hardness[0] / hardness[1];
			f = ground[src][1] + f2;
		}
		deltaSedimentGround[src] = f;
	} else {
		// depositing sediment
		deltaSedimentGround[src] = Kd * delta;
	}
	*/
}

template<bool safe>
inline void Grid::ErosionAndDepositionUpdate(int x, int y) {
	HydroPure::ErosionAndDepositionUpdate(x, y);
	/*
	int src = At<false>(x, y);
	float ds = deltaSedimentGround[src];
	ground[src].AddGeneral(-ds);
	sediment[src] += ds;
	deltaSedimentGround[src] = 0;
	*/
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	HydroPure::SedimentTransportation(x, y);
	/*
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	FOR_EACH_DIR_SAFE_COND(neighs[DIR],
		{
			const float sum = SumFlux(neighs[DIR]);
			if (sum > 0.0f) {
				const float neighSed = sediment[neighs[DIR]];
				const float incomingFlux = flux[neighs[DIR]].fluxArray[R_DIR];
				const float dV = dt * neighSed * incomingFlux / sum;
				deltaSedimentGround[src] += dV;
			}
		});
	*/
}

template<bool safe>
inline void Grid::SedimentTransportationUpdate(int x, int y) {
	HydroPure::SedimentTransportationUpdate(x, y);
	/*
	int src = At<false>(x, y);
	float f = SumFlux(src);
	if (f > 0.0f) {
		sediment[src] *= 1.0f - dt;
	}
	sediment[src] += deltaSedimentGround[src];
	deltaSedimentGround[src] = 0;
	*/
}

template<bool safe>
inline void Grid::ThermalErosionCalculation(int x, int y) {
	HydroPure::ThermalErosionCalculation(x, y);
	/*
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	NEIGHBOURS_CORNERS(neighCorners, x, y);
	deltaSedimentGround[src] = 0;
	float _h11 = ground[src][1];
	float _h0 = ground[src][0];
	float _h1 = _h0 + _h11;
	float t1 = l * tangentOfAngleOfRecluse[1];
	float t0 = l * tangentOfAngleOfRecluse[0];
	constexpr float halfSqrt2 = 0.7071067811865475244f;
	float delta = 0.0f;
	
	float n0;
	float n11;
	float n1;
	
	auto func = [&](){
		float h0 = _h0;
		float h1 = _h1;
		float h11 = _h11;

		bool swapped = h1 < n1;
		if (swapped) {
			std::swap(n0, h0);
			std::swap(n1, h1);
			std::swap(n11, h11);
		}

		float hh0 = n0 + t0;
		float hh1 = n1 + t1;
		float d = 0.0f;
		if (hh0 < hh1) {
			hh0 = hh1;
		}

		if (h1 > hh1) {
			d = std::min(h1 - hh1, h11);
			if (h1 > hh0) {
				d = std::max(d, h1 - hh0);
			}
			delta = delta + (swapped ? d : -d);
		}
	};
	
	FOR_EACH_DIR_SAFE_COND(neighs[DIR],
		{
			n0 = ground[neighs[DIR]][0];
			n11 = ground[neighs[DIR]][1];
			n1 = n0 + n11;
			func();
		});
	if (true) {
		t1 *= 2.0f * halfSqrt2;
		t0 *= 2.0f * halfSqrt2;
		FOR_EACH_DIR_SAFE_COND(neighCorners[DIR],
			{
				n0 = ground[neighCorners[DIR]][0];
				n11 = ground[neighCorners[DIR]][1];
				n1 = n0 + n11;
				func();
			});
	}
	deltaSedimentGround[src] = delta * (1.0f/8.0f) * 0.5f * dt;
	*/
}

template<bool safe>
inline void Grid::ThermalErosionUpdate(int x, int y) {
	HydroPure::ThermalErosionUpdate(x, y);
	/*
	int src = At<false>(x, y);
	ground[src].AddGeneral(deltaSedimentGround[src]);
	*/
}

inline float Grid::EvaporationRate(int x, int y) {
	return 0.04; // Make it dependent on temperature in place (x,y)
}

template<bool safe>
inline void Grid::Evaporation(int x, int y) {
	HydroPure::Evaporation(x, y);
	/*
	int src = At<false>(x, y);
	water[src] *= (1 - EvaporationRate(x, y)*dt);
	*/
}

template<bool safe>
inline void Grid::Smooth(int x, int y) {
	HydroPure::Smooth(x, y);
	/*
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	constexpr float mult = 128;
	float sum = ground[src].Total() * (mult - 4);
	FOR_EACH_DIR(sum += neighs[DIR] ? ground[neighs[DIR]].Total() : ground[src].Total());
	deltaSedimentGround[src] = sum / mult;
	*/
}

template<bool safe>
inline void Grid::SmoothUpdate(int x, int y) {
	HydroPure::SmoothUpdate(x, y);
	/*
	int src = At<false>(x, y);
	float dh = deltaSedimentGround[src] - ground[src].Total();
	ground[src][0] += dh;
	*/
}

template<bool safe>
inline void Grid::ClearDelta(int x, int y) {
	HydroPure::ClearDelta(x, y);
	/*
	int src = At<false>(x, y);
	deltaSedimentGround[src] = 0;
	*/
}

constexpr int XDXDX = 4;

static std::vector<std::thread> threads;
static std::atomic<int> jobId, jobsDone, jobsTotal;
static std::function<void(int X)> jobWorker;

template<int _BORDER, bool PARALLEL, typename T1, typename T2>
inline void Grid::ForEachSafeBorders(T1 &&funcSafe, T2 &&funcUnsafe)
{
	constexpr int BORDER = 1; // _BORDER;
	static std::function<void()> singleIteration =
		[]()
	{
		int X = jobId.fetch_add(1);
		if (X < jobsTotal.load()) {
			jobWorker(X);
			jobsDone++;
		}
	};
	{
		int threadsCount =
			std::max(((int)std::thread::hardware_concurrency()) - 2, 0) / 2;
		while (threads.size() < threadsCount) {
			threads.push_back(std::thread(
						[](){
						while(true) {
							if (jobId.load() < jobsTotal.load()) {
								singleIteration();
							} else {
								std::this_thread::sleep_for(std::chrono::milliseconds(1));
							}
						}
						}));
		}
	}
	const bool parallel = PARALLEL && !threads.empty();
	
	if (parallel) {
		jobWorker = [&](int X){
			X = X * XDXDX + BORDER;
			int xend = std::min(X+XDXDX, width-BORDER);
			int yend = height-1;
			for(int y=1; y<yend; ++y) {
				for (int x=X; x<xend; ++x) {
					funcUnsafe(x, y);
				}
			}
		};
		jobsDone = 0;
		jobId = 0;
		jobsTotal = (width-BORDER*2 + XDXDX - 1) / XDXDX;
	}
	
	if (!parallel) {
		for(int X=BORDER; X<width-BORDER; X+=XDXDX) {
			for(int y=1; y<height-1; ++y) {
				for (int x=X; x<X+XDXDX && x<width-BORDER; ++x) {
					funcUnsafe(x, y);
				}
			}
		}
	}
	
	for(int i=BORDER; i<width-BORDER; ++i) {
		for(int j=0; j<BORDER && j<height/2; ++j) {
			funcSafe(i, j);
			funcSafe(i, height-1-j);
		}
	}
	for(int i=0; i<height; ++i) {
		for(int j=0; j<BORDER && j<width/2; ++j) {
			funcSafe(j, i);
			funcSafe(width-1-j, i);
		}
	}
	
	if (parallel) {
		while (jobsDone.load() < jobsTotal.load()) {
			singleIteration();
		}
	}

	jobsTotal = 0;
	jobsDone = 0;
	jobId = 0;
}

#define FOR_EACH_SAFE_BORDERS(BORDER, PARALLEL, FUNC) \
	ForEachSafeBorders<BORDER, PARALLEL>([this](int x, int y){this->FUNC<true>(x, y);}, [this](int x, int y){this->FUNC<false>(x, y);});

void Grid::FullCycle() {
	HydroPure::water = water;
	HydroPure::ground = ground;
	HydroPure::sediment = sediment;
	HydroPure::deltaSedimentGround = deltaSedimentGround; 
	HydroPure::velocity = velocity;
	HydroPure::flux = flux;

	HydroPure::hardness[0] = hardness[0];
	HydroPure::hardness[1] = hardness[1];
	HydroPure::tangentOfAngleOfRecluse[0] = tangentOfAngleOfRecluse[0];
	HydroPure::tangentOfAngleOfRecluse[1] = tangentOfAngleOfRecluse[1];
	HydroPure::width = width;
	HydroPure::height = height;
	HydroPure::dt = dt;
	HydroPure::A = A;
	HydroPure::g = g;
	HydroPure::l = l;
	HydroPure::Kd = Kd;
	HydroPure::Kc = Kc;
	HydroPure::minimumSedimentCapacity = minimumSedimentCapacity;
	
	++iter;
	constexpr bool parallel = false;
	if (useWater) {
		FOR_EACH_SAFE_BORDERS(1, parallel, CalcOutflux);
		FOR_EACH_SAFE_BORDERS(1, parallel, UpdateWaterLevelAndVelocity);
		FOR_EACH_SAFE_BORDERS(1, parallel, ErosionAndDepositionCalculation);
		FOR_EACH_SAFE_BORDERS(0, parallel, ErosionAndDepositionUpdate);
// 		FOR_EACH_SAFE_BORDERS(0, false, ClearDelta);
		FOR_EACH_SAFE_BORDERS(1, parallel, SedimentTransportation);
		FOR_EACH_SAFE_BORDERS(0, parallel, SedimentTransportationUpdate);
	}
	if (useThermalErosion) {
		FOR_EACH_SAFE_BORDERS(1, parallel, ThermalErosionCalculation);
		FOR_EACH_SAFE_BORDERS(0, parallel, ThermalErosionUpdate);
	}
	if (useWater) {
		FOR_EACH_SAFE_BORDERS(0, parallel, Evaporation);
	}
	if (useSmoothing && iter % 47 == 0) {
		// TODO: replace with selectional smoothing, to smooth only where slope
		// changes very rapidly
		FOR_EACH_SAFE_BORDERS(1, parallel, Smooth);
		FOR_EACH_SAFE_BORDERS(0, parallel, SmoothUpdate);
	}
}

#undef SAFE_COND_GRID
