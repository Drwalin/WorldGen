#pragma once

#include <chrono>
#include <cmath>
#include <cassert>

#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <functional>

#ifndef HYDRO_EROSION_INCL_HPP
#define HYDRO_EROSION_INCL_HPP

#include "HydroErosion.hpp"

#include "HydroErosionMacros.hpp"

template<bool safe>
inline Tile* Grid::At(int x, int y) const {
	if constexpr(safe) {
// 		if(x < 0)            return At<true>(width-1, y);//nullptr;
// 		else if(x >= width)  return At<true>(0, y);//nullptr;
// 		if(y < 0)            return At<true>(x, height-1);//nullptr;
// 		else if(y >= height) return At<true>(x, 0);//nullptr;
		if(x < 0)            return nullptr;
		else if(x >= width)  return nullptr;
		if(y < 0)            return nullptr;
		else if(y >= height) return nullptr;
	}
	return (Tile*)(tiles + (x*height + y));
}

template<bool safe, int dir>
inline Tile* Grid::Neighbour(int x, int y) const {
	if constexpr(dir == 0) return At<safe>(x-1, y);
	if constexpr(dir == 1) return At<safe>(x, y+1);
	if constexpr(dir == 2) return At<safe>(x+1, y);
	if constexpr(dir == 3) return At<safe>(x, y-1);
	static_assert(dir >= 0 && dir <= 3);
}

template<bool safe>
inline Tile* Grid::Neighbour(int x, int y, int dir) const {
	switch(dir) {
		case 0: return Neighbour<safe, 0>(x, y);
		case 1: return Neighbour<safe, 1>(x, y);
		case 2: return Neighbour<safe, 2>(x, y);
		case 3: return Neighbour<safe, 3>(x, y);
	}
	assert(dir >= 0 && dir <= 3);
	// TODO: should not happen
}

template<bool safe, int dir>
inline float Grid::CalcFluxInDirection(Tile& src, Tile* neigh) const {
	if(neigh == nullptr)
		return 0.0f;
	Tile& dst = *neigh;
	float dh = (src.b - dst.b) + (src.d - dst.d);
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
	src.d += (dt/(l*l)) * fs;
	if (src.d < 0.0f) {
		src.d = 0;
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
	float ds = src.deltaSedimentGround * 0.25;
	src.ground -= ds;
	src.sediment += ds;
	src.deltaSedimentGround = 0;
	
	if (src.sediment < 0.0f) {
		printf("SDUP(%.16f)\n", src.sediment);
	}
}

static inline float SumFlux(const Tile &t) {
	return t.f.L + t.f.B + t.f.R + t.f.T;
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	Tile& src = *At<false>(x, y);
	if constexpr (false) {
		float sx, sy;
		sx = (x - src.vx*dt)/l;
		sy = (y - src.vy*dt)/l;
		
		int SX = sx;
		int SY = sy;
		
		if(SX < 0)
			return;
		if(SX+1 >= width)
			return;
		if(SY < 0)
			return;
		if(SY+1 >= height)
			return;
		
		Tile& a00 = *At<false>(SX, SY);
		Tile& a10 = *At<false>(SX+1, SY);
		Tile& a11 = *At<false>(SX+1, SY+1);
		Tile& a01 = *At<false>(SX, SY+1);
		
		float rx = sx - SX;
		float ry = sy - SY;
		
		float ds00 = a00.s * (1-rx) * (1-ry);
		float ds10 = a10.s * (rx) * (1-ry);
		float ds01 = a01.s * (1-rx) * (ry);
		float ds11 = a11.s * (rx) * (ry);
		
		a00.s -= ds00;
		a10.s -= ds10;
		a01.s -= ds01;
		a11.s -= ds11;
		
		src.s += ds00 + ds10 + ds01 + ds11;
	} else {
		/*
		src.deltaSedimentGround = 0;
		
		Tile *tiles[3][3];
		for (int X = 0; X <= 2; ++X) {
			for (int Y = 0; Y <= 2; ++Y) {
				tiles[X][Y] = At<safe>(x+X-1, y+Y-1);
			}
		}
		
		for (int _X = 0; _X <= 2; ++_X) {
			for (int _Y = 0; _Y <= 2; ++_Y) {
				if (_X == 1 && _Y == 1) {
					continue;
				}
				Tile *t = tiles[_X][_Y];
				int X = _X+x-1;
				int Y = _Y+y-1;
				float vx = t->vx;
				float vy = t->vy;
				float avx = vx < 0.0f ? -vx : vx;
				float avy = vy < 0.0f ? -vy : vy;
				
				if (vx * (_X-1.0f) >= 0.0f) {
					continue;
				}
				if (vy * (_Y-1.0f) >= 0.0f) {
					continue;
				}
				
				if (avx > avy) {
					
				} else if (avx == avy) {
					if (abs(_X-1) == abs(_Y-1)) {
						continue;
					}
					src.deltaSedimentGround += t->sediment;
					
					
				} else {
					
				}
				
				
				
				
				
				
				if (sum > 0.0f) {
					const float neighSed = neighs[DIR]->s;
					const float incomingFlux = neighs[DIR]->fluxArray[R_DIR];
					src.deltaSedimentGround += neighSed * incomingFlux / sum;
				}
				
			}
		}
		*/
		
		
		
		NEIGHBOURS(neighs, x, y);
		src.deltaSedimentGround = 0;
		float sum;
		FOR_EACH_DIR_SAFE_COND(neighs[DIR],
			{
				sum = SumFlux(*(neighs[DIR]));
				if (sum > 0.0f) {
					const float neighSed = neighs[DIR]->s;
					const float incomingFlux = neighs[DIR]->fluxArray[R_DIR];
					src.deltaSedimentGround += neighSed * incomingFlux / sum;
				}
			});
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
	
	
	float f = SumFlux(src);
	if (f > 0.0f) {
		src.s = src.deltaSedimentGround;
	} else {
		src.s += src.deltaSedimentGround;
	}
	src.deltaSedimentGround = 0;
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
	constexpr float mult = 512;
	float sum = src.ground * (mult - 4);
	FOR_EACH_DIR(sum += neighs[DIR] ? neighs[DIR]->ground : src.ground);
	src.deltaSedimentGround = sum / mult;
}

template<bool safe>
inline void Grid::SmoothUpdate(int x, int y) {
	Tile& src = *At<false>(x, y);
	src.ground = src.deltaSedimentGround;
}

constexpr int XDXDX = 4;

static std::vector<std::thread> threads;
static std::atomic<int> jobId, jobsDone, jobsTotal;
static std::function<void(int X)> jobWorker;

template<int BORDER, typename T1, typename T2>
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
	
// 	for(int X=BORDER; X<width-BORDER; X+=XDXDX) {
// 		for(int y=1; y<height-1; ++y) {
// 			for (int x=X; x<X+XDXDX && x<width-BORDER; ++x) {
// 				funcUnsafe(x, y);
// 			}
// 		}
// 	}
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
	
	while (jobsDone.load() < jobsTotal.load()) {
		singleIteration();
	}
	
	jobsTotal = 0;
	jobsDone = 0;
	jobId = 0;
}

#define FOR_EACH_SAFE_BORDERS(BORDER, FUNC) \
	ForEachSafeBorders<BORDER>([this](int x, int y){this->FUNC<true>(x, y);}, [this](int x, int y){this->FUNC<false>(x, y);});
/*
	for(int X=BORDER; X<width-BORDER; X+=XDXDX) { \
		for(int y=1; y<height-1; ++y) { \
			for (int x=X; x<X+XDXDX && x<width-BORDER; ++x) { \
				FUNC<false>(x, y); \
			} \
		} \
	} \
	for(int i=BORDER; i<width-BORDER; ++i) { \
		for(int j=0; j<BORDER && j<height/2; ++j) { \
			FUNC<true>(i, j); \
			FUNC<true>(i, height-1-j); \
		} \
	} \
	for(int i=0; i<height; ++i) { \
		for(int j=0; j<BORDER && j<width/2; ++j) { \
			FUNC<true>(j, i); \
			FUNC<true>(width-1-j, i); \
		} \
	}
	*/

inline void Grid::FullCycle() {
	++iter;
 	FOR_EACH_SAFE_BORDERS(1, CalcOutflux);
 	FOR_EACH_SAFE_BORDERS(1, UpdateWaterLevelAndVelocity);
 	FOR_EACH_SAFE_BORDERS(1, ErosionAndDepositionCalculation);
 	FOR_EACH_SAFE_BORDERS(0, ErosionAndDepositionUpdate);
 	FOR_EACH_SAFE_BORDERS(1, SedimentTransportation);
 	FOR_EACH_SAFE_BORDERS(0, SedimentTransportationUpdate);
 	FOR_EACH_SAFE_BORDERS(0, Evaporation);
	if (false && iter % 16 == 0) {
		FOR_EACH_SAFE_BORDERS(1, Smooth);
		FOR_EACH_SAFE_BORDERS(0, SmoothUpdate);
	}
}

#undef SAFE_COND_GRID

#endif

