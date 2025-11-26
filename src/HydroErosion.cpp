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
inline float Grid::CalcFluxInDirection(Tile& src, Tile* neigh) const {
	if(neigh == nullptr)
		return 0.0f;
	Tile& dst = *neigh;
	float dh = (src.ground - dst.ground) + (src.sediment - dst.sediment) + (src.water - dst.water);
	float f = src.fluxArray[dir] + dt * A * g * dh / l;
	if(f < 0.0f)
		f = 0.0f;
	return f;
}

inline void Grid::LimitFlux(Tile& src) {
	const float outflux = src.f.L + src.f.B + src.f.R + src.f.T;
	const float water = src.d * l*l;
	if(outflux <= 0.001) {
		FOR_EACH_DIR({src.fluxArray[DIR] = 0.0f;});
		return;
	}
	float K = water / (outflux * dt);
	if (K > 1.0f) {
		K = 1.0f;
	}
	if (K >= 0.0f) {
		FOR_EACH_DIR({src.fluxArray[DIR] *= K;});
	}
}

template<bool safe>
void Grid::CalcOutflux(int x, int y) {
	Tile& src = *At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], (src.fluxArray[DIR] = CalcFluxInDirection<safe, DIR>(src, neighs[DIR])));
	LimitFlux(src);
}

template<bool safe>
void Grid::UpdateWaterLevel(Tile& src, Tile** neighs) {
	float fs = 0;
	FOR_EACH_DIR_SAFE_COND(neighs[DIR], fs += neighs[DIR]->fluxArray[R_DIR] - src.fluxArray[DIR]);
	src.water += (dt/(l*l)) * fs;
	if (src.water < 0.0f) {
		src.water = 0;
	}
}

template<bool safe>
void Grid::UpdateWaterLevelAndVelocity(int x, int y) {
	Tile& src = *At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	float water_level = src.d;
	UpdateWaterLevel<safe>(src, neighs);
	water_level = (water_level + src.d) * 0.5f;
	if(water_level < 0.001) {
		src.vx = 0;
		src.vy = 0;
		return;
	}
	float dWx = 0.0f, dWy = 0.0f;
	SAFE_COND_GRID(neighs[0], dWx -= neighs[0]->fluxArray[2] - src.fluxArray[0]);
	SAFE_COND_GRID(neighs[1], dWy -= neighs[1]->fluxArray[3] - src.fluxArray[1]);
	SAFE_COND_GRID(neighs[2], dWx += neighs[2]->fluxArray[0] - src.fluxArray[2]);
	SAFE_COND_GRID(neighs[3], dWy += neighs[3]->fluxArray[1] - src.fluxArray[3]);
	src.vx = dWx / (water_level * l);
	src.vy = dWy / (water_level * l);
	
// 	if ( x < 9 && y < 9 && x > 5 && y > 5) {
// 		if ( x == 6 && y == 6) {
// 			printf("\n\n");
// 		}
// 		printf("dWX = %f     dWY = %f     wl*l = %f\n", dWx, dWy, water_level*l);
// 	}
}

template<bool safe>
inline float Grid::SinusLocalTiltAngle(Tile& t, int x, int y) {
	float xl, yl;
	NEIGHBOURS(neighs, x, y);
	if constexpr (safe) {
		xl = yl = 2.0f * l;
		if(neighs[0] == nullptr) { neighs[0] = &t; xl = l; }
		if(neighs[1] == nullptr) { neighs[1] = &t; yl = l; }
		if(neighs[2] == nullptr) { neighs[2] = &t; xl = l; }
		if(neighs[3] == nullptr) { neighs[3] = &t; yl = l; }
	} else {
		xl = yl = l;
	}
	
	const float dhdx = (neighs[0]->b - neighs[2]->b) / xl;
	const float dhdy = (neighs[1]->b - neighs[3]->b) / yl;
	
	const float s = dhdx*dhdx + dhdy*dhdy;
	
	return sqrt(s) / sqrt(1.0f + s);
}

template<bool safe>
inline void Grid::ErosionAndDepositionCalculation(int x, int y) {
	Tile& src = *At<false>(x, y);
	const float sinusLocalTiltAngle = SinusLocalTiltAngle<safe>(src, x, y);
	const float v = sqrt(src.vx * src.vx + src.vy * src.vy);
	float capacity = Kc * sinusLocalTiltAngle * v;
	if(capacity < minimumSedimentCapacity)
		capacity = minimumSedimentCapacity;
	if (capacity > 1.0f)
		capacity = 1.0f;
	
	const float delta = capacity - src.sediment;
	
	if(capacity > src.sediment) {
		// picking up sediment
		src.deltaSedimentGround = src.Ks * delta;
		if (src.sediment + src.deltaSedimentGround < 0.0f && src.sediment >= 0.0f) {
			printf(" PICK cap = %f   sed = %f    srcDelta = %f     delta = %f     result = %f\n", capacity,
					src.sediment, src.deltaSedimentGround, delta, src.sediment + src.deltaSedimentGround);
		}
	} else {
		// depositing sediment
		src.deltaSedimentGround = Kd * delta;
		if (src.sediment + src.deltaSedimentGround < 0.0f && src.sediment >= 0.0f) {
			printf(" DEP cap = %f   sed = %f    srcDelta = %f     delta = %f     result = %f\n", capacity,
					src.sediment, src.deltaSedimentGround, delta, src.sediment + src.deltaSedimentGround);
		}
	}
}

template<bool safe>
inline void Grid::ErosionAndDepositionUpdate(int x, int y) {
	Tile& src = *At<false>(x, y);
	float ds = src.deltaSedimentGround * 0.9f;
	src.ground -= ds;
	src.sediment += ds;
	src.deltaSedimentGround = 0;
	
// 	if (src.sediment < 0.0f) {
// 		printf("SDUP(%.16f)\n", src.sediment);
// 	}
}

static inline float SumFlux(const Tile &t) {
	return t.f.L + t.f.B + t.f.R + t.f.T;
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	Tile& src = *At<false>(x, y);
	switch(3) {
	case 0: {
		const float sx = (x - src.vx*dt)/l;
		const float sy = (y - src.vy*dt)/l;
		
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
		
		Tile empty;
		empty.sediment = 0;
		Tile *a00 = At<safe>(SX, SY);
		Tile *a10 = At<safe>(SX+1, SY);
		Tile *a11 = At<safe>(SX+1, SY+1);
		Tile *a01 = At<safe>(SX, SY+1);
		if (a00 == nullptr) a00 = &empty;
		if (a10 == nullptr) a10 = &empty;
		if (a01 == nullptr) a01 = &empty;
		if (a11 == nullptr) a11 = &empty;
		
		const float rx = sx - SX;
		const float ry = sy - SY;
		
		const float ds00 = 0.25f * a00->s * (1.0f-rx) * (1.0f-ry);
		const float ds10 = 0.25f * a10->s * (rx) * (1.0f-ry);
		const float ds01 = 0.25f * a01->s * (1.0f-rx) * (ry);
		const float ds11 = 0.25f * a11->s * (rx) * (ry);
		
		const static auto F = +[](float v)->bool{ return v >= 0.0f && v <= 1.0f; };
		assert(F(rx));
		assert(F(1.0f-rx));
		assert(F(ry));
		assert(F(1.0f-ry));
		
		a00->s -= ds00;
		a10->s -= ds10;
		a01->s -= ds01;
		a11->s -= ds11;
		
		src.s += ds00 + ds10 + ds01 + ds11;
	} break;
	case 1: {
		NEIGHBOURS(neighs, x, y);
		FOR_EACH_DIR_SAFE_COND(neighs[DIR],
			{
				const float sum = SumFlux(*(neighs[DIR])) * 16.0f;
				if (sum > 0.0f) {
					const float neighSed = neighs[DIR]->s;
					const float incomingFlux = neighs[DIR]->fluxArray[R_DIR];
					const float dV = neighSed * incomingFlux / sum;
					src.deltaSedimentGround += dV;
					neighs[DIR]->deltaSedimentGround -= dV;
				}
			});
	} break;
	case 2: {
		float vx = src.vx;
		float vy = src.vy;
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
		
		Tile *nx = At<safe>(x+dvx, y);
		Tile *ny = At<safe>(x, y+dvy);
		Tile *nxy = At<safe>(x+dvx, y+dvy);
		
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
		
		const float sediment = src.sediment;
		const float toDiag = sediment * diagv * factor;
		const float toX = sediment * fvx * factor;
		const float toY = sediment * fvy * factor;
		
		const float toTransfer = toDiag + toX + toY;
		
		src.deltaSedimentGround -= toTransfer;
		if (nx) nx->deltaSedimentGround += toX;
		if (ny) ny->deltaSedimentGround += toY;
		if (nxy) nxy->deltaSedimentGround += toDiag;
		
		assert (toTransfer <= sediment + 0.0000001);
		
	} break;
	case 3: {
		const float sx = (x - src.vx*dt)/l;
		const float sy = (y - src.vy*dt)/l;
		
		const int SX = floor(sx);
		const int SY = floor(sy);
		
		Tile empty;
		empty.sediment = 0;
		Tile *a00 = At<safe>(SX, SY);
		Tile *a10 = At<safe>(SX+1, SY);
		Tile *a11 = At<safe>(SX+1, SY+1);
		Tile *a01 = At<safe>(SX, SY+1);
		if (a00 == nullptr) a00 = &empty;
		if (a10 == nullptr) a10 = &empty;
		if (a01 == nullptr) a01 = &empty;
		if (a11 == nullptr) a11 = &empty;
		
		const float rx = sx - SX;
		const float ry = sy - SY;
		
		const float ds00 = a00->s * (1.0f-rx) * (1.0f-ry);
		const float ds10 = a10->s * (rx) * (1.0f-ry);
		const float ds01 = a01->s * (1.0f-rx) * (ry);
		const float ds11 = a11->s * (rx) * (ry);
		
		const static auto F = +[](float v)->bool{ return v >= 0.0f && v <= 1.0f; };
		assert(F(rx));
		assert(F(1.0f-rx));
		assert(F(ry));
		assert(F(1.0f-ry));
		
		a00->s -= ds00;
		a10->s -= ds10;
		a01->s -= ds01;
		a11->s -= ds11;
		
		src.deltaSedimentGround = (ds00 + ds10 + ds01 + ds11);
	} break;
	}
}

template<bool safe>
inline void Grid::SedimentTransportationUpdate(int x, int y) {
	Tile& src = *At<false>(x, y);
	
	
	for (int i=0; i<4; ++i) {
		float f = src.fluxArray[i];
		if (f < 1000000 && f > -1000000) {
		} else {
			printf("FLUX(%i,%.16f)", i, f);
		}
	}
	
	
// 	float f = SumFlux(src);
// 	if (f > 0.0f) {
// 		src.s = src.deltaSedimentGround;
// 	} else {
		src.s += src.deltaSedimentGround;
// 	}
	src.deltaSedimentGround = 0;
	
	
	if (src.sediment < 0.0f) {
		printf("SDUP2(%.16f)\n", src.sediment);
	}
}

inline float Grid::EvaporationRate(int x, int y) {
	return 0.01; // Make it dependent on temperature in place (x,y)
}

template<bool safe>
inline void Grid::Evaporation(int x, int y) {
	Tile& src = *At<false>(x, y);
	src.d *= (1 - EvaporationRate(x, y)*dt);
}

template<bool safe>
inline void Grid::Smooth(int x, int y) {
	Tile& src = *At<false>(x, y);
	NEIGHBOURS(neighs, x, y);
	constexpr float mult = 128;
	float sum = src.ground * (mult - 4);
	FOR_EACH_DIR(sum += neighs[DIR] ? neighs[DIR]->ground : src.ground);
	src.deltaSedimentGround = sum / mult;
}

template<bool safe>
inline void Grid::SmoothUpdate(int x, int y) {
	Tile& src = *At<false>(x, y);
	src.ground = src.deltaSedimentGround;
}

template<bool safe>
inline void Grid::ClearDelta(int x, int y) {
	Tile& src = *At<false>(x, y);
	src.deltaSedimentGround = 0;
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
			std::max(((int)std::thread::hardware_concurrency()) - 2, 0);
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
			for(int y=1; y<height-1; ++y) {
				for (int x=X; x<X+XDXDX && x<width-BORDER; ++x) {
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
 	FOR_EACH_SAFE_BORDERS(0, false, ClearDelta);
 	FOR_EACH_SAFE_BORDERS(12, false, SedimentTransportation);
 	FOR_EACH_SAFE_BORDERS(0, false, SedimentTransportationUpdate);
 	FOR_EACH_SAFE_BORDERS(0, true, Evaporation);
	if (true && iter % 1 == 0) {
		FOR_EACH_SAFE_BORDERS(1, true, Smooth);
		FOR_EACH_SAFE_BORDERS(0, true, SmoothUpdate);
	}
}

#undef SAFE_COND_GRID
