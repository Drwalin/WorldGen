
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
		x = x%width;
		y = y%height;
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

template TileId HydroErosion::At<false>(int x, int y) const;
template TileId HydroErosion::At<true>(int x, int y) const;

template<bool safe, int dir>
TileId HydroErosion::Neighbour(int x, int y) const {
	if constexpr(dir == 0) return At<safe>(x-1, y);
	if constexpr(dir == 1) return At<safe>(x, y-1);
	if constexpr(dir == 2) return At<safe>(x+1, y);
	if constexpr(dir == 3) return At<safe>(x, y+1);
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

template<bool safe>
void HydroErosion::GetNeighbours(int x, int y, TileId* neighs) const {
	neighs[0] = Neighbour<safe, 0>(x, y);
	neighs[1] = Neighbour<safe, 1>(x, y);
	neighs[2] = Neighbour<safe, 2>(x, y);
	neighs[3] = Neighbour<safe, 3>(x, y);
}

template<bool safe>
void HydroErosion::GetNeighboursOrSelf(int x, int y, TileId self, TileId* neighs) const {
	neighs[0] = Neighbour<safe, 0>(x, y);
	neighs[1] = Neighbour<safe, 1>(x, y);
	neighs[2] = Neighbour<safe, 2>(x, y);
	neighs[3] = Neighbour<safe, 3>(x, y);
	for(int i=0; i<4; ++i) {
		if(neighs[i] < 0) {
			neighs[i] = self;
		}
	}
}

template<bool safe>
void HydroErosion::GetNeighboursOrSelf(int x, int y, TileId* neighs) const {
	GetNeighboursOrSelf<safe>(x, y, At<safe>(x, y), neighs);
}




template<bool safe>
void HydroErosion::CalcFluxInDirection(TileId src, TileId neigh, int dir) {
	const float dh = b[src] + d[src] - b[neigh] - d[neigh];
	fluxArray[src].v[dir] = std::max<float>(
			0.f,
			fluxArray[src].v[dir] + dt * A * g * dh / l
			);
}

void HydroErosion::LimitFlux(TileId src) {
	float water = d[src] * l * l;
	float outflux = (f[src].L + f[src].B + f[src].R + f[src].T) * dt;
	if(outflux < 0.0001) {
		f[src].L = 0;
		f[src].B = 0;
		f[src].R = 0;
		f[src].T = 0;
		return;
	}
	if(water < outflux) {
		const float K = std::min(1.f, water / outflux);
		f[src].L *= K;
		f[src].B *= K;
		f[src].R *= K;
		f[src].T *= K;
	}
}

// 3.2.1
template<bool safe>
void HydroErosion::CalcOutflux(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4];
	GetNeighbours<safe>(x, y, neighs);
	for(int i=0; i<4; ++i)
		CalcFluxInDirection<safe>(src, neighs[i], i);
	LimitFlux(src);
}





// 3.2.2
template<bool safe>
void HydroErosion::UpdateWaterLevelAndVelocity(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4];
	GetNeighbours<safe>(x, y, neighs);
	
	const float DV = (
		  fluxArray[neighs[0]].v[2] - fluxArray[src].v[0]
		+ fluxArray[neighs[1]].v[3] - fluxArray[src].v[1]
		+ fluxArray[neighs[2]].v[0] - fluxArray[src].v[2]
		+ fluxArray[neighs[3]].v[1] - fluxArray[src].v[3]
		) * dt;
	
	float d_ = d[src];
	
	d[src] += DV / (l*l);
	
	d_ = (d[src] + d_)*0.5f;
	if(d_ < 0.001) {
		vx[src] = 0;
		vy[src] = 0;
		return;
	}
	
	
	
	
	const float dWx = fluxArray[neighs[0]].v[2] - fluxArray[src].v[0]
		- fluxArray[neighs[2]].v[0] + fluxArray[src].v[2];
	const float dWy = -fluxArray[neighs[1]].v[3] + fluxArray[src].v[1]
		+ fluxArray[neighs[3]].v[1] - fluxArray[src].v[3];
	
	vx[src] = dWx / (l * d_) / 2;
	vy[src] = dWy / (l * d_) / 2;
}





template<bool safe>
float HydroErosion::SinusLocalTiltAngle(TileId src, int x, int y) const {
	TileId neighs[4];
	GetNeighboursOrSelf<safe>(x, y, neighs);
	float h[4];
	for(int i=0; i<4; ++i) {
		h[i] = ground[neighs[i]];
	}
	float dbdx = (h[2]-h[0]) / (2.f * l);
	float dbdy = (h[1]-h[3]) / (2.f * l);
	float dxy = dbdx*dbdx + dbdy*dbdy;
	return sqrt(dxy) / sqrt(1 + dxy);
}

// 3.3
template<bool safe>
void HydroErosion::ErosionAndDeposition(int x, int y) {
	TileId src = At<false>(x, y);
	const float sin_tilt = std::max<float>(SinusLocalTiltAngle<safe>(src, x, y), 0.0);
	float C = Kc * sin_tilt * sqrt(vx[src]*vx[src] + vy[src]*vy[src]);
	
	if(C > sediment[src]) { // deposit sediment to ground
		const float amount = Ks * (C - sediment[src]) * dt;
		ground[src] -= amount;
		sediment[src] += amount;
	} else if(C < sediment[src]) { // pick up soil from ground into deposit
		const float amount = Kd * (sediment[src] - C) * dt;
		ground[src] += amount;
		sediment[src] -= amount;
	}
}





// 3.4
template<bool safe>
void HydroErosion::SedimentTransportation(int x, int y) {
	return ;
	
	TileId dst = At<false>(x, y);
	float sx, sy;
	sx = x - vx[dst]*dt/l;
	sy = y - vy[dst]*dt/l;
	
	int SX = sx;
	int SY = sy;
	
	TileId a00 = At<false>(SX, SY);
	TileId a10 = At<false>(SX+1, SY);
	TileId a11 = At<false>(SX+1, SY+1);
	TileId a01 = At<false>(SX, SY+1);
	
	const float rx = sx - (float)SX;
	const float ry = sy - (float)SY;
	
	const float ds00 = s[a00];
	const float ds10 = s[a10];
	const float ds11 = s[a11];
	const float ds01 = s[a01];
	
	const float s_top = rx * ds10 + (1-rx) * ds00;
	const float s_bot = rx * ds11 + (1-rx) * ds01;
	
	newSediment[dst] = ry * s_top + (1-ry) * s_bot;
	
	
	
// 	TileId dst = At<false>(x, y);
// 	float sx, sy;
// 	sx = (x*l + vx[dst]*dt)/l;
// 	sy = (y*l + vy[dst]*dt)/l;
// 	
// 	int SX = sx;
// 	int SY = sy;
// 	
// 	TileId a00 = At<safe>(SX, SY);
// 	TileId a10 = At<safe>(SX+1, SY);
// 	TileId a11 = At<safe>(SX+1, SY+1);
// 	TileId a01 = At<safe>(SX, SY+1);
// 	
// 	const float rx = sx - (float)SX;
// 	const float ry = sy - (float)SY;
// 	
// 	const float sed = s[dst];
// 	s[dst] = 0;
// 	
// 	float ds00 = sed * (1-rx) * (1-ry);
// 	float ds10 = sed * (rx) * (1-ry);
// 	float ds01 = sed * (1-rx) * (ry);
// 	float ds11 = sed * (rx) * (ry);
// 	
// 	newSediment[a00] += ds00;
// 	newSediment[a10] += ds10;
// 	newSediment[a01] += ds01;
// 	newSediment[a11] += ds11;
}

template<bool safe>
void HydroErosion::SedimentTransportationUpdateWhatsLeft(int x, int y) {
	TileId src = At<safe>(x, y);
	s[src] = newSediment[src];
	newSediment[src] = 0;
}





float HydroErosion::EvaporationRate(int x, int y) {
	return 0.1;
}

// 3.5
template<bool safe>
void HydroErosion::Evaporation(int x, int y) {
	TileId src = At<false>(x, y);
	d[src] = std::max(0.f, d[src] * (1.f - 0.1f*dt));
}





#define FOR_EACH_SAFE_BORDERS(BORDER, FUNC) \
	for(int x=BORDER; x<width-BORDER; ++x) { \
		for(int y=BORDER; y<height-BORDER; ++y) { \
			FUNC<false>(x, y); \
		} \
	} \
	for(int i=0; i<width; ++i) { \
		for(int j=0; j<BORDER; ++j) { \
			FUNC<true>(i, j); \
			FUNC<true>(i, height-1-j); \
		} \
	} \
	for(int i=BORDER; i<height-BORDER; ++i) { \
		for(int j=0; j<BORDER; ++j) { \
			FUNC<true>(j, i); \
			FUNC<true>(width-1-j, i); \
		} \
	}

// to be executed after water increase
void HydroErosion::FullCycle() {
	for(int i=0; i<width*height; ++i) {
		newSediment[i] = 0;
	}
	struct three {
		float water;
		float sediment;
		float ground;
	};
	auto SumWater = [this](){
		float sumWater = 0;
		for(int i=0; i<width*height; ++i) {
			sumWater += water[i];
		}
		float sumGround = 0;
		for(int i=0; i<width*height; ++i) {
			sumGround += ground[i];
		}
		float sumSed = 0;
		for(int i=0; i<width*height; ++i) {
			sumSed += s[i];
		}
		return three{sumWater/(width*height), sumSed/(width*height), sumGround/(width*height)};
	};
	
	auto PrintWater = [&]() {
		auto w = SumWater();
		printf("water: %f,  sediment: %f,  ground: %f   all dirt: %f\n", w.water, w.sediment, w.ground, w.sediment+w.ground);
	};
		
	
	printf("\n\n");
	PrintWater();
 	FOR_EACH_SAFE_BORDERS((wrap?0:1), CalcOutflux);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS((wrap?0:1), UpdateWaterLevelAndVelocity);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS((wrap?0:4), ErosionAndDeposition);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS(0, SedimentTransportation);
 	FOR_EACH_SAFE_BORDERS(0, SedimentTransportationUpdateWhatsLeft);
	PrintWater();
	
// 	std::swap(newSediment, sediment);
	
 	FOR_EACH_SAFE_BORDERS(0, Evaporation);
}





