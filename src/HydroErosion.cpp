#include <cmath>
#include <cassert>

#include <chrono>
#include <vector>
#include <thread>
#include <atomic>
#include <functional>

#include "../include/worldgen/HydroErosion.hpp"
#include "../include/worldgen/HydroErosionMacros.hpp"

std::atomic<float> f;

int a = (f += 1.31, 0);

template<bool safe, int dir>
inline float Grid::CalcFluxInDirection(int src, int neigh) const {
	if(neigh == 0)
		return 0.0f;
	int dst = neigh;
	float dh = (ground[src] - ground[dst]) + (sediment[src] - sediment[dst]) + (water[src] - water[dst]);
	float f = flux[src].fluxArray[dir] + dt * A * g * dh / l;
	if(f < 0.0f)
		f = 0.0f;
	return f;
}

inline void Grid::LimitFlux(int src) {
	const float outflux = flux[src].f.L + flux[src].f.B + flux[src].f.R + flux[src].f.T;
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

template<bool safe>
void Grid::CalcOutflux(int x, int y) {
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], (flux[src].fluxArray[DIR] = CalcFluxInDirection<safe, DIR>(src, neighs[DIR])));
	LimitFlux(src);
}

template<bool safe>
void Grid::UpdateWaterLevel(int src, int* neighs) {
	float fs = 0;
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], fs += flux[neighs[DIR]].fluxArray[R_DIR] - flux[src].fluxArray[DIR]);
	water[src] += (dt/(l*l)) * fs;
	if (water[src] < 0.0f) {
		water[src] = 0;
	}
}

template<bool safe>
void Grid::UpdateWaterLevelAndVelocity(int x, int y) {
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
	
// 	if ( x < 9 && y < 9 && x > 5 && y > 5) {
// 		if ( x == 6 && y == 6) {
// 			printf("\n\n");
// 		}
// 		printf("dWX = %f     dWY = %f     wl*l = %f\n", dWx, dWy, water_level*l);
// 	}
}

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
	
	const float dhdx = (ground[neighs[0]] - ground[neighs[2]]) / xl;
	const float dhdy = (ground[neighs[1]] - ground[neighs[3]]) / yl;
	const float s = dhdx*dhdx + dhdy*dhdy;
	return sqrt(s) / sqrt(1.0f + s);
}

template<bool safe>
inline void Grid::ErosionAndDepositionCalculation(int x, int y) {
	int src = At<false>(x, y);
	const float sinusLocalTiltAngle = SinusLocalTiltAngle<safe>(src, x, y);
	const float v = sqrt(velocity[src].x * velocity[src].x + velocity[src].y * velocity[src].y);
	float capacity = Kc * sinusLocalTiltAngle * v;
	if(capacity < minimumSedimentCapacity)
		capacity = minimumSedimentCapacity;
	if (capacity > 1.0f)
		capacity = 1.0f;
	
	const float delta = capacity - sediment[src];
	
	if(capacity > sediment[src]) {
		// picking up sediment
		deltaSedimentGround[src] = hardness[src] * delta;
		if (sediment[src] + deltaSedimentGround[src] < 0.0f && sediment[src] >= 0.0f) {
			printf(" PICK cap = %f   sed = %f    srcDelta = %f     delta = %f     result = %f\n", capacity,
					sediment[src], deltaSedimentGround[src], delta, sediment[src] + deltaSedimentGround[src]);
		}
	} else {
		// depositing sediment
		deltaSedimentGround[src] = Kd * delta;
		if (sediment[src] + deltaSedimentGround[src] < 0.0f && sediment[src] >= 0.0f) {
			printf(" DEP cap = %f   sed = %f    srcDelta = %f     delta = %f     result = %f\n", capacity,
					sediment[src], deltaSedimentGround[src], delta, sediment[src] + deltaSedimentGround[src]);
		}
	}
}

template<bool safe>
inline void Grid::ErosionAndDepositionUpdate(int x, int y) {
	int src = At<false>(x, y);
	float ds = deltaSedimentGround[src] * dt;
	ground[src] -= ds;
	sediment[src] += ds;
	deltaSedimentGround[src] = 0;
	
// 	if (src.sediment < 0.0f) {
// 		printf("SDUP(%.16f)\n", src.sediment);
// 	}
}

inline float Grid::SumFlux(int t) {
	return flux[t].f.L + flux[t].f.B + flux[t].f.R + flux[t].f.T;
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	int src = At<false>(x, y);
	switch(5) {
	case 0: {
		const float sx = (x - velocity[src].x*dt)/l;
		const float sy = (y - velocity[src].y*dt)/l;
		
		const int SX = floor(sx);
		const int SY = floor(sy);
		
// 		if(SX < 0)
// 			return;
// 		if(SX+1 >= width)
// 			return;
// 		if(SY < 0)
// 			return;
// 		if(SY+1 >= height)
// 			return;
		
		int a00 = At<safe>(SX, SY);
		int a10 = At<safe>(SX+1, SY);
		int a11 = At<safe>(SX+1, SY+1);
		int a01 = At<safe>(SX, SY+1);
		
		const float rx = sx - SX;
		const float ry = sy - SY;
		
		const float ds00 = 0.25f * (a00==0 ? 0.0f : sediment[a00]) * (1.0f-rx) * (1.0f-ry);
		const float ds10 = 0.25f * (a10==0 ? 0.0f : sediment[a10]) * (rx) * (1.0f-ry);
		const float ds01 = 0.25f * (a01==0 ? 0.0f : sediment[a01]) * (1.0f-rx) * (ry);
		const float ds11 = 0.25f * (a11==0 ? 0.0f : sediment[a11]) * (rx) * (ry);
		
		const static auto F = +[](float v)->bool{ return v >= 0.0f && v <= 1.0f; };
		assert(F(rx));
		assert(F(1.0f-rx));
		assert(F(ry));
		assert(F(1.0f-ry));
		
		if (a00) sediment[a00] -= ds00;
		if (a10) sediment[a10] -= ds10;
		if (a01) sediment[a01] -= ds01;
		if (a11) sediment[a11] -= ds11;
		
		sediment[src] += ds00 + ds10 + ds01 + ds11;
	} break;
	case 1: {
		NEIGHBOURS(neighs, x, y);
		FOR_EACH_DIR_SAFE_COND(neighs[DIR],
			{
				const float sum = SumFlux(neighs[DIR]) * 16.0f;
				if (sum > 0.0f) {
					const float neighSed = sediment[neighs[DIR]];
					const float incomingFlux = flux[neighs[DIR]].fluxArray[R_DIR];
					const float dV = neighSed * incomingFlux / sum;
					deltaSedimentGround[src] += dV;
					deltaSedimentGround[neighs[DIR]] -= dV;
				}
			});
	} break;
	case 2: {
		float vx = velocity[src].x;
		float vy = velocity[src].y;
		float avx = vx < 0.0f ? -vx : vx;
		float avy = vy < 0.0f ? -vy : vy;
		
		const float maxv = avx > avy ? avx : avy;
		const float minv = avx < avy ? avx : avy;
		
		if (maxv > 0.0001f) {
		} else {
			break;
		}
		
		const int dvx = vx > 0.0f ? 1 : (vx < 0.0f ? -1 : 0);
		const int dvy = vy > 0.0f ? 1 : (vy < 0.0f ? -1 : 0);
		
		int nx = At<safe>(x+dvx, y);
		int ny = At<safe>(x, y+dvy);
		int nxy = At<safe>(x+dvx, y+dvy);
		
		const float _diagv = nxy ? minv : 0.0f;
		float fvx = nx ? avx - _diagv : 0.0f;
		float fvy = ny ? avy - _diagv : 0.0f;
		float diagv = _diagv * 1.4142135623730950488; // sqrt(2)
		
		float total = diagv + fvx + fvy;
		if (total > 0.95f) {
			float div = 1.0f / (total * 1.05f);
// 			diagv /= total;
// 			fvx /= total;
// 			fvy /= total;
			diagv *= div;
			fvx *= div;
			fvy *= div;
			total = 1.0f;
			if (nx && ny && nxy) {
// 				assert(fabs(total - (diagv + fvx + fvy)) < 0.00001);
			}
		}
		
		if (nx && ny && nxy) {
// 			assert(fabs(total - (diagv + fvx + fvy)) < 0.00001);
		} else {
			total = diagv + fvx + fvy;
		}
		
// 		assert(total <= 1.0f);
		
		constexpr float factor = 0.125f;
		
		const float sediment = this->sediment[src];
		const float toDiag = sediment * diagv * factor;
		const float toX = sediment * fvx * factor;
		const float toY = sediment * fvy * factor;
		
		const float toTransfer = toDiag + toX + toY;
		
		deltaSedimentGround[src] -= toTransfer;
		if (nx) deltaSedimentGround[nx] += toX;
		if (ny) deltaSedimentGround[ny] += toY;
		if (nxy) deltaSedimentGround[nxy] += toDiag;
		
		assert (toTransfer <= sediment + 0.0000001);
		
	} break;
	case 3: {
		const float sx = (x - velocity[src].x*dt)/l;
		const float sy = (y - velocity[src].y*dt)/l;
		
		const int SX = floor(sx);
		const int SY = floor(sy);
		
		int a00 = At<safe>(SX, SY);
		int a10 = At<safe>(SX+1, SY);
		int a11 = At<safe>(SX+1, SY+1);
		int a01 = At<safe>(SX, SY+1);
		
		const float rx = sx - SX;
		const float ry = sy - SY;
		
		const float ds00 = a00==0 ? 0.0f : sediment[a00] * (1.0f-rx) * (1.0f-ry);
		const float ds10 = a10==0 ? 0.0f : sediment[a10] * (rx) * (1.0f-ry);
		const float ds01 = a01==0 ? 0.0f : sediment[a01] * (1.0f-rx) * (ry);
		const float ds11 = a11==0 ? 0.0f : sediment[a11] * (rx) * (ry);
		
		const static auto F = +[](float v)->bool{ return v >= 0.0f && v <= 1.0f; };
		assert(F(rx));
		assert(F(1.0f-rx));
		assert(F(ry));
		assert(F(1.0f-ry));
		
		if (a00) sediment[a00] -= ds00;
		if (a10) sediment[a10] -= ds10;
		if (a01) sediment[a01] -= ds01;
		if (a11) sediment[a11] -= ds11;
		
		deltaSedimentGround[src] = (ds00 + ds10 + ds01 + ds11);
	} break;
	case 4: {
		const float sx = (x - velocity[src].x*dt)/l;
		const float sy = (y - velocity[src].y*dt)/l;
		const float vel = sqrt(velocity[src].x*velocity[src].x + velocity[src].y*velocity[src].y);
		const float capFactor = (vel - 1.0f) * 0.001f * dt;
		if (capFactor < 0.0f) {
			break;
		}
		const float capacity = capFactor > 0.1f ? 0.1f : capFactor;
		
		const int SX = floor(sx);
		const int SY = floor(sy);
		
		int a00 = At<safe>(SX, SY);
		int a10 = At<safe>(SX+1, SY);
		int a11 = At<safe>(SX+1, SY+1);
		int a01 = At<safe>(SX, SY+1);
		
		const float rx = sx - SX;
		const float ry = sy - SY;
		
		const float ds00 = a00 == 0 ? 0.0f : capacity * (1.0f-rx) * (1.0f-ry);
		const float ds10 = a10 == 0 ? 0.0f : capacity * (rx) * (1.0f-ry);
		const float ds01 = a01 == 0 ? 0.0f : capacity * (1.0f-rx) * (ry);
		const float ds11 = a11 == 0 ? 0.0f : capacity * (rx) * (ry);
		
		if (a00) { ground[a00] -= ds00; ground[src] += ds00; }
		if (a10) { ground[a10] -= ds10; ground[src] += ds10; }
		if (a01) { ground[a01] -= ds01; ground[src] += ds01; }
		if (a11) { ground[a11] -= ds11; ground[src] += ds11; }
	} break;
	case 5: {
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
	} break;
	}
}

template<bool safe>
inline void Grid::SedimentTransportationUpdate(int x, int y) {
	int src = At<false>(x, y);
	
	
	for (int i=0; i<4; ++i) {
		float f = flux[src].fluxArray[i];
		if (f < 1000000 && f > -1000000) {
		} else {
			printf("FLUX(%i,%.16f)", i, f);
		}
	}
	
	float f = SumFlux(src);
	if (f > 0.0f) {
		sediment[src] *= 1.0f - dt;
	}
	sediment[src] += deltaSedimentGround[src];
	
	/*
	float f = SumFlux(src);
	if (f > 0.0f) {
		src.s = src.deltaSedimentGround;
	} else {
		src.sediment += src.deltaSedimentGround;
	}
	*/
	deltaSedimentGround[src] = 0;
	
	
// 	if (src.sediment < 0.0f) {
// 		printf("SDUP2(%.16f)\n", src.sediment);
// 	}
}

inline float Grid::EvaporationRate(int x, int y) {
	return 0.01; // Make it dependent on temperature in place (x,y)
}

template<bool safe>
inline void Grid::Evaporation(int x, int y) {
	int src = At<false>(x, y);
	water[src] *= (1 - EvaporationRate(x, y)*dt);
}

template<bool safe>
inline void Grid::Smooth(int x, int y) {
	int src = At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	constexpr float mult = 128;
	float sum = ground[src] * (mult - 4);
	FOR_EACH_DIR(sum += neighs[DIR] ? ground[neighs[DIR]] : ground[src]);
	deltaSedimentGround[src] = sum / mult;
}

template<bool safe>
inline void Grid::SmoothUpdate(int x, int y) {
	int src = At<false>(x, y);
	ground[src] = deltaSedimentGround[src];
}

template<bool safe>
inline void Grid::ClearDelta(int x, int y) {
	int src = At<false>(x, y);
	deltaSedimentGround[src] = 0;
}

constexpr int XDXDX = 4;

static std::vector<std::thread> threads;
static std::atomic<int> jobId, jobsDone, jobsTotal;
static std::function<void(int X)> jobWorker;

template<int BORDER, bool PARALLEL, typename T1, typename T2>
inline void Grid::ForEachSafeBorders(T1 &&funcSafe, T2 &&funcUnsafe)
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
	++iter;
	FOR_EACH_SAFE_BORDERS(1, true, CalcOutflux);
	FOR_EACH_SAFE_BORDERS(1, true, UpdateWaterLevelAndVelocity);
	FOR_EACH_SAFE_BORDERS(1, true, ErosionAndDepositionCalculation);
	FOR_EACH_SAFE_BORDERS(0, true, ErosionAndDepositionUpdate);
// 	FOR_EACH_SAFE_BORDERS(0, false, ClearDelta);
	FOR_EACH_SAFE_BORDERS(1, true, SedimentTransportation);
	FOR_EACH_SAFE_BORDERS(0, true, SedimentTransportationUpdate);
	FOR_EACH_SAFE_BORDERS(0, true, Evaporation);
	if (false && iter % 1 == 0) {
		FOR_EACH_SAFE_BORDERS(1, true, Smooth);
		FOR_EACH_SAFE_BORDERS(0, true, SmoothUpdate);
	}
}

#undef SAFE_COND_GRID
