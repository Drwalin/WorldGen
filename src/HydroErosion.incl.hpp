
#pragma once

#ifndef HYDRO_EROSION_INCL_HPP
#define HYDRO_EROSION_INCL_HPP

#include "HydroErosion.hpp"


#define SAFE_COND_GRID(COND, EXPR) \
	if constexpr (safe) { \
		if(COND) { \
			EXPR; \
		} \
	} else { \
		EXPR; \
	}

template<bool safe>
inline Tile* Grid::At(int x, int y) const {
	if constexpr(safe) {
		if(x < 0)            return At<true>(width-1, y);//NULL;
		else if(x >= width)  return At<true>(0, y);//NULL;
		if(y < 0)            return At<true>(x, height-1);//NULL;
		else if(y >= height) return At<true>(x, 0);//NULL;
	}
	return (Tile*)(tiles + (x*height + y));
}

template<bool safe, int dir>
inline Tile* Grid::Neighbour(int x, int y) const {
	if constexpr(dir == 0) return At<safe>(x-1, y);
	if constexpr(dir == 1) return At<safe>(x, y+1);
	if constexpr(dir == 2) return At<safe>(x+1, y);
	if constexpr(dir == 3) return At<safe>(x, y-1);
}

template<bool safe>
inline Tile* Grid::Neighbour(int x, int y, int dir) const {
	switch(dir) {
		case 0: return Neighbour<safe, 0>(x, y);
		case 1: return Neighbour<safe, 1>(x, y);
		case 2: return Neighbour<safe, 2>(x, y);
		case 3: return Neighbour<safe, 3>(x, y);
	}
	// TODO: should not happen
}

template<bool safe, int dir>
inline float Grid::CalcFluxInDirection(Tile& src, Tile& neigh) const {
	Tile* dstp = &neigh;
	if(dstp == NULL)
		return 0.0f;
	Tile& dst = *dstp;
	float dh = src.b + src.d - dst.b - dst.d;
	float f = src.fluxArray[dir] +  dt * A * g * dh / l;
	if(f < 0)
		f = 0;
	return f;
}

inline void Grid::LimitFlux(Tile& src) {
	float sum = src.f.L + src.f.B + src.f.R + src.f.T;
	float water = src.d * l*l;
	float outflux = sum * dt;
	if(outflux <= water+0.000001)
		return;
	float K = water / outflux;
	src.f.L *= K;
	src.f.B *= K;
	src.f.R *= K;
	src.f.T *= K;
}

template<bool safe>
inline void Grid::CalcOutflux(int x, int y) {
	Tile& src = *At<false>(x, y);
	Tile* neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	SAFE_COND_GRID(true, (src.fluxArray[0] = CalcFluxInDirection<safe, 0>(src, *(neighs[0]))));
	SAFE_COND_GRID(true, (src.fluxArray[1] = CalcFluxInDirection<safe, 1>(src, *(neighs[1]))));
	SAFE_COND_GRID(true, (src.fluxArray[2] = CalcFluxInDirection<safe, 2>(src, *(neighs[2]))));
	SAFE_COND_GRID(true, (src.fluxArray[3] = CalcFluxInDirection<safe, 3>(src, *(neighs[3]))));
	LimitFlux(src);
}

template<bool safe>
inline void Grid::UpdateWaterLevel(Tile& src, Tile** neighs) {
	float fs = 0;
	SAFE_COND_GRID(neighs[0], fs += neighs[0]->f.R - src.f.L);
	SAFE_COND_GRID(neighs[1], fs += neighs[1]->f.L - src.f.R);
	SAFE_COND_GRID(neighs[2], fs += neighs[2]->f.T - src.f.B);
	SAFE_COND_GRID(neighs[3], fs += neighs[3]->f.B - src.f.T);
	src.d += dt/l/l * fs;
}

template<bool safe>
inline void Grid::UpdateWaterLevelAndVelocity(int x, int y) {
	Tile& src = *At<false>(x, y);
	Tile* neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	float water_level = src.d;
	UpdateWaterLevel<safe>(src, neighs);
	water_level = (water_level + src.d) * 0.5f;
	if(water_level < 0.000001) {
		src.vx = 0;
		src.vy = 0;
		return;
	}
	float dWx = 0.0f, dWy = 0.0f;
	SAFE_COND_GRID(neighs[0], dWx += neighs[0]->f.R - src.f.L);
	SAFE_COND_GRID(neighs[2], dWx -= neighs[2]->f.L - src.f.R);
	SAFE_COND_GRID(neighs[1], dWy += neighs[1]->f.T - src.f.B);
	SAFE_COND_GRID(neighs[3], dWy -= neighs[3]->f.B - src.f.T);
	src.vx = dWx / (water_level * l);
	src.vy = dWy / (water_level * l);
}

template<bool safe>
inline float Grid::SinusLocalTiltAngle(Tile& t, int x, int y) {
	float xl, yl;
	Tile* neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	if constexpr (safe) {
		xl = yl = 2.0f * l;
		if(neighs[0] == NULL) { neighs[0] = &t; xl = l; }
		if(neighs[1] == NULL) { neighs[1] = &t; yl = l; }
		if(neighs[2] == NULL) { neighs[2] = &t; xl = l; }
		if(neighs[3] == NULL) { neighs[3] = &t; yl = l; }
	} else {
		xl = yl = l;
	}
	
	float dhdx = (neighs[0]->b - neighs[2]->b) / xl;
	float dhdy = (neighs[1]->b - neighs[3]->b) / yl;
	
	float dhdx2 = dhdx*dhdx;
	float dhdy2 = dhdy*dhdy;
	
	float s = dhdx2 + dhdy2;
	
	return sqrt(s) / sqrt(1.0f + s);
}

template<bool safe>
inline void Grid::ErosionAndDeposition(int x, int y) {
	Tile& src = *At<false>(x, y);
	float sinusLocalTiltAngle = SinusLocalTiltAngle<safe>(src, x, y);
	float C = Kc * sinusLocalTiltAngle * sqrt(src.vx*src.vx + src.vy*src.vy);
	if(C < minimumSedimentCapacity)
		C = minimumSedimentCapacity;
	if(C > src.s) {
		float d = src.Ks*(C - src.s);
		// Here src.b should be updated in next step for whole ground instead in another phase to prevent broblems with calculating tilt angle before updating sediment
		src.b -= d;
		src.s += d;
	} else {
		float d = Kd*(src.s - C);
		// Here src.b should be updated in next step
		src.b += d;
		src.s -= d;
	}
}

template<bool safe>
inline void Grid::SedimentTransportation(int x, int y) {
	Tile& src = *At<false>(x, y);
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
}

inline float Grid::EvaporationRate(int x, int y) {
	return 0.01; // Make it dependent on temperature in place (x,y)
}

template<bool safe>
inline void Grid::Evaporation(int x, int y) {
	Tile& src = *At<false>(x, y);
	src.d *= (1 - EvaporationRate(x, y)*dt);
}

#define FOR_EACH_SAFE_BORDERS(BORDER, FUNC) \
	for(int x=BORDER; x<width-BORDER; ++x) { \
		for(int y=1; y<height-1; ++y) { \
			FUNC<false>(x, y); \
		} \
	} \
	for(int i=0; i<width; ++i) { \
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

inline void Grid::FullCycle() {
 	FOR_EACH_SAFE_BORDERS(1, CalcOutflux);
 	FOR_EACH_SAFE_BORDERS(1, UpdateWaterLevelAndVelocity);
 	FOR_EACH_SAFE_BORDERS(1, ErosionAndDeposition);
 	FOR_EACH_SAFE_BORDERS(16, SedimentTransportation);
 	FOR_EACH_SAFE_BORDERS(1, Evaporation);
}

#undef FOR_EACH_SAFE_BORDERS

#undef SAFE_COND_GRID

#endif

