
#include <algorithm>

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
TileId HydroErosion::At(int x, int y) const {
	if(wrap) {
		return ((x+width)%width) + ((y+height)%height) * width;
	} else {
		if constexpr(safe) {
			if(x < 0)            return -1;
			else if(x >= width)  return -1;
			if(y < 0)            return -1;
			else if(y >= height) return -1;
		}
		return x + y * width;
	}
}

template<bool safe, int dir>
TileId HydroErosion::Neighbour(int x, int y) const {
	if constexpr(dir == 0) return At<safe>(x-1, y);
	if constexpr(dir == 1) return At<safe>(x, y+1);
	if constexpr(dir == 2) return At<safe>(x+1, y);
	if constexpr(dir == 3) return At<safe>(x, y-1);
}

template<bool safe>
TileId HydroErosion::Neighbour(int x, int y, int dir) const {
	switch(dir) {
		case 0: return Neighbour<safe, 0>(x, y);
		case 1: return Neighbour<safe, 1>(x, y);
		case 2: return Neighbour<safe, 2>(x, y);
		case 3: return Neighbour<safe, 3>(x, y);
	}
}





template<bool safe, int dir>
float HydroErosion::CalcFluxInDirection(TileId src, TileId neigh) const {
	const float f_old = fluxArray[src].v[dir];
	const float dh = b[src] + d[src] - b[neigh] - d[neigh];
	return std::max<float>(
			0,
			f_old + dt * A * g * dh / l
			);
}

void HydroErosion::LimitFlux(TileId src) {
	float sum = f[src].L + f[src].B + f[src].R + f[src].T;
	float water = d[src] * l*l;
	float outflux = sum * dt;
	if(outflux <= water+0.000001 || outflux <= 0.000001)
		return;
	float K = water / outflux;
	f[src].L *= K;
	f[src].B *= K;
	f[src].R *= K;
	f[src].T *= K;
}

// 3.2.1
template<bool safe>
void HydroErosion::CalcOutflux(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	SAFE_COND_GRID(neighs[0]>=0, (fluxArray[src].v[0] = CalcFluxInDirection<safe, 0>(src, neighs[0])));
	SAFE_COND_GRID(neighs[1]>=0, (fluxArray[src].v[1] = CalcFluxInDirection<safe, 1>(src, neighs[1])));
	SAFE_COND_GRID(neighs[2]>=0, (fluxArray[src].v[2] = CalcFluxInDirection<safe, 2>(src, neighs[2])));
	SAFE_COND_GRID(neighs[3]>=0, (fluxArray[src].v[3] = CalcFluxInDirection<safe, 3>(src, neighs[3])));
	LimitFlux(src);
}





template<bool safe>
void HydroErosion::UpdateWaterLevel(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	float fs = 0;
	SAFE_COND_GRID(neighs[0]>=0, fs += f[neighs[0]].R - f[src].L);
	SAFE_COND_GRID(neighs[1]>=0, fs += f[neighs[1]].L - f[src].R);
	SAFE_COND_GRID(neighs[2]>=0, fs += f[neighs[2]].T - f[src].B);
	SAFE_COND_GRID(neighs[3]>=0, fs += f[neighs[3]].B - f[src].T);
	DV[src] = fs*dt/(l*l);
	d[src] += DV[src];
}

// 3.2.2
template<bool safe>
void HydroErosion::UpdateWaterVelocity(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	
	float dWx=0, dWy=0;
	
	
	
}





template<bool safe>
float HydroErosion::SinusLocalTiltAngle(TileId t, int x, int y) {
	return 0;
}

// 3.3
template<bool safe>
void HydroErosion::ErosionAndDeposition(int x, int y) {
}





// 3.4
template<bool safe>
void HydroErosion::SedimentTransportation(int x, int y) {
}





float HydroErosion::EvaporationRate(int x, int y) {
	return 0;
}

// 3.5
template<bool safe>
void HydroErosion::Evaporation(int x, int y) {
}





// to be executed after water increase
void HydroErosion::FullCycle() {
}





