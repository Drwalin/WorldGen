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

template<bool safe>
void Grid::CalcOutflux(int x, int y) {
	HydroPure::CalcOutflux(x, y);
}

template<bool safe>
void Grid::UpdateWaterLevelAndVelocity(int x, int y) {
	HydroPure::UpdateWaterLevelAndVelocity(x, y);
}

template<bool safe>
inline void Grid::ErosionAndDepositionCalculation(int x, int y) {
	HydroPure::ErosionAndDepositionCalculation(x, y);
}

template<bool safe>
inline void Grid::ErosionAndDepositionUpdate(int x, int y) {
	HydroPure::ErosionAndDepositionUpdate(x, y);
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	HydroPure::SedimentTransportation(x, y);
}

template<bool safe>
inline void Grid::SedimentTransportationUpdate(int x, int y) {
	HydroPure::SedimentTransportationUpdate(x, y);
}

template<bool safe>
inline void Grid::ThermalErosionCalculation(int x, int y) {
	HydroPure::ThermalErosionCalculation(x, y);
}

template<bool safe>
inline void Grid::ThermalErosionUpdate(int x, int y) {
	HydroPure::ThermalErosionUpdate(x, y);
}

template<bool safe>
inline void Grid::Evaporation(int x, int y) {
	HydroPure::Evaporation(x, y);
}

template<bool safe>
inline void Grid::Smooth(int x, int y) {
	HydroPure::Smooth(x, y);
}

template<bool safe>
inline void Grid::SmoothUpdate(int x, int y) {
	HydroPure::SmoothUpdate(x, y);
}

template<bool safe>
inline void Grid::ClearDelta(int x, int y) {
	HydroPure::ClearDelta(x, y);
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
