#include <cassert>

#include <chrono>
#include <vector>
#include <thread>
#include <atomic>
#include <functional>

#include "../include/worldgen/HydroErosion.hpp"
#include "HydroErosionPure.h"

void Grid::CallHydroErosion()
{
	if (useGpu) {
		gpu.CallUpdateRainAndRiver();
		gpu.CallCalcOutFlux();
		gpu.CallUpdateWaterLevelAndVelocity();
		gpu.CallErosionAndDepositionCalculation();
		gpu.CallErosionAndDepositionUpdate();
		gpu.CallSedimentTransportation();
		gpu.CallSedimentTransportationUpdate();
// 		gpu.vboSediment->Copy(gpu.vboTemp1, 0, 0, elements * 4);
	} else {
		ForEachSafeBorders(HydroPure::CalcOutFlux);
		ForEachSafeBorders(HydroPure::UpdateWaterLevelAndVelocity);
		ForEachSafeBorders(HydroPure::ErosionAndDepositionCalculation);
		ForEachSafeBorders(HydroPure::ErosionAndDepositionUpdate);
		ForEachSafeBorders(HydroPure::SedimentTransportation);
		ForEachSafeBorders(HydroPure::SedimentTransportationUpdate);
// 		memcpy(sediment, temp1, elements * 4);
	}
}

void Grid::CallThermalErosion()
{
	if (useGpu) {
		gpu.CallThermalErosionCalculation();
// 		gpu.vboGround->Copy(gpu.vboVelocity, 0, 0, elements * 8);
		gpu.CallThermalErosionUpdate();
	} else {
		ForEachSafeBorders(HydroPure::ThermalErosionCalculation);
		ForEachSafeBorders(HydroPure::ThermalErosionUpdate);
// 		memcpy(ground, velocity, elements * 8);
	}
}
void Grid::CallEvaporation() {
	if (useGpu) {
		gpu.CallEvaporation();
		gpu.CallEvaporationUpdate();
// 		gpu.vboWater->Copy(gpu.vboTemp1, 0, 0, elements * 4);
	} else {
		ForEachSafeBorders(HydroPure::Evaporation);
		ForEachSafeBorders(HydroPure::EvaporationUpdate);
// 		memcpy(water, temp1, elements * 4);
	}
}
void Grid::CallSmoothing()
{
	// TODO: replace with selectional smoothing, to smooth only where slope
	// changes very rapidly
	if (useGpu) {
		gpu.CallSmooth();
		gpu.CallSmoothUpdate();
	} else {
		ForEachSafeBorders(HydroPure::Smooth);
		ForEachSafeBorders(HydroPure::SmoothUpdate);
	}
}

void Grid::Init(int w, int h, bool useGpu)
{
	width = w;
	height = h;
	if (w & (ALIGNEMENT-1) || h & (ALIGNEMENT-1)) {
		printf("Grid width and height need to be multiple of %i\n", ALIGNEMENT);
		fflush(stdout);
		exit(1);
	}
	
	elements = width * height + 1;
	elementsStorage = elements + ALIGNEMENT;
	this->useGpu = useGpu;
	
	ground = new GroundLayers[elementsStorage] + OFFSET;
	velocity = new Velocity[elementsStorage] + OFFSET;
	flux = new Flux[elementsStorage] + OFFSET;
	water_sediment_temp = new glm::vec4[elementsStorage];
	
	water = ((float*)water_sediment_temp) + OFFSET;
	sediment = water + elementsStorage;
	temp1 = sediment + elementsStorage;
	
	for (int i = 0; i < elements; ++i) {
		water[i] = 0.0f;
		sediment[i] = 0.0f;
		temp1[i] = 0.0f;
		flux[i].f[0] = 0;
		flux[i].f[1] = 0;
		flux[i].f[2] = 0;
		flux[i].f[3] = 0;
		ground[i].layers[0] = 0;
		ground[i].layers[1] = 0;
	}
	if (useGpu) {
		gpu.Init(width, height, this);
	}
}

Grid::Grid()
{
	parallelThreads = std::max(((int)std::thread::hardware_concurrency()) - 2, 0) / 2;
	width = height = 0;
	dt = 0.03;
	crossSectionalAreaOfPipe = .6;
	gravity = 9.81;
	tileDimensionSize = 1;

	depositionConstant = 0.03;
	sedimentCapacityConstant = 0.01;
	minimumSedimentCapacity = 0.01;
}

Grid::~Grid()
{
	if (ground) {
		delete[] (ground - OFFSET);
	}
	if (velocity) {
		delete[] (velocity - OFFSET);
	}
	if (flux) {
		delete[] (flux - OFFSET);
	}
	if (water_sediment_temp) {
		delete water_sediment_temp;
	}
}

constexpr int XDXDX = 4;

static std::vector<std::thread> threads;
static std::atomic<int> jobId, jobsDone, jobsTotal;
static std::function<void(int X)> jobWorker;

template <typename TFunc> inline void Grid::ForEachSafeBorders(TFunc &&func)
{
	parallelThreads = std::clamp<int>(parallelThreads, 1, std::thread::hardware_concurrency());
	static std::function<void()> singleIteration = []() {
		int X = jobId.fetch_add(1);
		if (X < jobsTotal.load()) {
			jobWorker(X);
			jobsDone++;
		}
	};
	{
		while (threads.size() < parallelThreads-1) {
			threads.push_back(std::thread([]() {
				while (true) {
					if (jobId.load() < jobsTotal.load()) {
						singleIteration();
					} else {
						std::this_thread::sleep_for(
							std::chrono::milliseconds(1));
					}
				}
			}));
		}
	}
	const bool parallel = this->parallel && !threads.empty();

	if (parallel) {
		jobWorker = [&](int X) {
			X = X * XDXDX;
			int xend = std::min(X + XDXDX, width);
			int yend = height;
			for (int y = 0; y < yend; ++y) {
				for (int x = X; x < xend; ++x) {
					func(x, y);
				}
			}
		};
		jobsDone = 0;
		jobId = 0;
		jobsTotal = (width + XDXDX - 1) / XDXDX;

		if (parallel) {
			while (jobsDone.load() < jobsTotal.load()) {
				singleIteration();
			}
		}
	} else {
		for (int X = 0; X < width; X += XDXDX) {
			for (int y = 0; y < height; ++y) {
				for (int x = X; x < X + XDXDX && x < width; ++x) {
					func(x, y);
				}
			}
		}
	}

	jobsTotal = 0;
	jobsDone = 0;
	jobId = 0;
}

void Grid::FullCycle()
{
	HydroPure::water = water;
	HydroPure::ground = ground;
	HydroPure::sediment = sediment;
	HydroPure::temp1 = temp1;
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
	
	HydroPure::iteration = iteration;

	if (useWater) {
		CallHydroErosion();
	}
	if (useThermalErosion) {
		CallThermalErosion();
	}
	if (useWater) {
		CallEvaporation();
	}
	if (useSmoothing && iteration % 47 == 0) {
		// TODO: replace with selectional smoothing, to smooth only where slope
		// changes very rapidly
		CallSmoothing();
	}
	
	++iteration;
}

void Grid::UpdateHeightsTexture(gl::Texture *tex)
{
	gpu.UpdateHeightsTexture(tex);
}

#undef SAFE_COND_GRID
