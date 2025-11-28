// #include <cmath>
#include <cassert>

#include <chrono>
#include <vector>
#include <thread>
#include <atomic>
#include <functional>

#include "../include/worldgen/HydroErosion.hpp"
#include "HydroErosionPure.hpp"

constexpr int XDXDX = 4;

static std::vector<std::thread> threads;
static std::atomic<int> jobId, jobsDone, jobsTotal;
static std::function<void(int X)> jobWorker;

template<typename TFunc>
inline void Grid::ForEachSafeBorders(TFunc &&func)
{
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
	const bool parallel = this->parallel && !threads.empty();
	
	if (parallel) {
		jobWorker = [&](int X){
			X = X * XDXDX;
			int xend = std::min(X+XDXDX, width);
			int yend = height;
			for(int y=0; y<yend; ++y) {
				for (int x=X; x<xend; ++x) {
					func(x, y);
				}
			}
		};
		jobsDone = 0;
		jobId = 0;
		jobsTotal = (width + XDXDX - 1) / XDXDX;
	}
	
	if (!parallel) {
		for(int X=0; X<width; X+=XDXDX) {
			for(int y=0; y<height; ++y) {
				for (int x=X; x<X+XDXDX && x<width; ++x) {
					func(x, y);
				}
			}
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
	if (useWater) {
		ForEachSafeBorders(HydroPure::CalcOutflux);
		ForEachSafeBorders(HydroPure::UpdateWaterLevelAndVelocity);
		ForEachSafeBorders(HydroPure::ErosionAndDepositionCalculation);
		ForEachSafeBorders(HydroPure::ErosionAndDepositionUpdate);
		ForEachSafeBorders(HydroPure::ClearDelta);
		ForEachSafeBorders(HydroPure::SedimentTransportation);
		ForEachSafeBorders(HydroPure::SedimentTransportationUpdate);
	}
	if (useThermalErosion) {
		ForEachSafeBorders(HydroPure::ThermalErosionCalculation);
		ForEachSafeBorders(HydroPure::ThermalErosionUpdate);
	}
	if (useWater) {
		ForEachSafeBorders(HydroPure::Evaporation);
	}
	if (useSmoothing && iter % 47 == 0) {
		// TODO: replace with selectional smoothing, to smooth only where slope
		// changes very rapidly
		ForEachSafeBorders(HydroPure::Smooth);
		ForEachSafeBorders(HydroPure::SmoothUpdate);
	}
}

#undef SAFE_COND_GRID
