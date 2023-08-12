
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
	const float dh = b[src] + s[src] + d[src] - b[neigh] - s[neigh] - d[neigh];
	return std::max<float>(
			0,
			f_old + dt * A * g * dh / l
			);
}

void HydroErosion::LimitFlux(TileId src) {
	float sum = f[src].L + f[src].B + f[src].R + f[src].T;
	float water = d[src] * A;
	float outflux = sum * dt;
	if(outflux < 0.0000001) {
		f[src].L = 0;
		f[src].B = 0;
		f[src].R = 0;
		f[src].T = 0;
		return;
	}
	if(outflux-0.0000001 < water) {
		return;
	}
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





// 3.2.2
template<bool safe>
void HydroErosion::UpdateWaterLevelAndVelocity(int x, int y) {
	TileId src = At<false>(x, y);
	TileId neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	
	float dWx=0, dWy=0;
	float lx=0, ly=0;
	SAFE_COND_GRID(neighs[0]>=0, (lx += l, dWx += fluxArray[neighs[0]].v[2] - fluxArray[src].v[0]));
	SAFE_COND_GRID(neighs[2]>=0, (lx += l, dWx += fluxArray[neighs[2]].v[0] - fluxArray[src].v[2]));
	SAFE_COND_GRID(neighs[1]>=0, (ly += l, dWy += fluxArray[neighs[1]].v[3] - fluxArray[src].v[1]));
	SAFE_COND_GRID(neighs[3]>=0, (ly += l, dWy += fluxArray[neighs[3]].v[1] - fluxArray[src].v[3]));
	const float fs = dWx + dWy;
	const float DV_A = fs*dt/(lx*ly);
	d[src] += DV_A;
	float d_ = d[src] - DV_A*0.5f;
	if(d_ < 0.000001) {
		vx[src] = 0;
		vy[src] = 0;
		return;
	}
	vx[src] = dWx / (lx * d_);
	vy[src] = dWy / (ly * d_);
}





template<bool safe>
float HydroErosion::SinusLocalTiltAngle(TileId src, int x, int y) const {
	float nx, ny, nz;
	float vax=0, vaz=0;
	float vby=0, vbz=0;
	const float vay=0, vbx=0;
	
	TileId neighs[4] = {
		Neighbour<safe, 0>(x, y),
		Neighbour<safe, 1>(x, y),
		Neighbour<safe, 2>(x, y),
		Neighbour<safe, 3>(x, y) };
	
	float xaz, xbz, yaz, ybz;
	xaz = xbz = yaz = ybz = d[src];
	
	SAFE_COND_GRID(neighs[0]>=0, (vax += l, xaz = b[neighs[0]]));
	SAFE_COND_GRID(neighs[2]>=0, (vax += l, xbz = b[neighs[2]]));
	SAFE_COND_GRID(neighs[1]>=0, (vby += l, yaz = b[neighs[1]]));
	SAFE_COND_GRID(neighs[3]>=0, (vby += l, ybz = b[neighs[3]]));
	
	vaz = xbz - xaz;
	vbz = ybz - yaz;
	
	nx = vay*vbz - vaz*vby;
	ny = vax*vbz - vaz*vbx;
	nz = vay*vbx - vax*vby;
	
	if(nz < 0) {
		nx = -nx;
		ny = -ny;
		nz = -nz;
	}
	
	const float nl = sqrt(nx*nx + ny*ny + nz*nz);
	
	nx /= nl;
	ny /= nl;
	nz /= nl;
	
	// cross product between unit normal and unit up
	
	const float tx = ny*1.f - nz*0.f;
	const float ty = nx*1.f - nz*0.f;
	const float tz = ny*0.f - nx*0.f;
	
	float ret = sqrt(tx*tx + ty*ty + tz*tz);
	return ret;
}

// 3.3
template<bool safe>
void HydroErosion::ErosionAndDeposition(int x, int y) {
	TileId src = At<false>(x, y);
	const float sin_tilt = std::max<float>(SinusLocalTiltAngle<safe>(src, x, y), 0.1f);
	const float C = Kc * sin_tilt * sqrt(vx[src]*vx[src] + vy[src]*vy[src]);
	
	if(C > s[src]) { // deposit sediment to ground
		const float amount = Ks * (C - s[src]);
		b[src] -= amount;
		s[src] += amount;
	} else { // pick up soil from ground into deposit
		const float amount = Kd * (s[src] - C);
		b[src] += amount;
		s[src] -= amount;
	}
}





// 3.4
template<bool safe>
void HydroErosion::SedimentTransportation(int x, int y) {
	TileId dst = At<false>(x, y);
	float sx, sy;
	sx = (x - vx[dst]*dt)/l;
	sy = (y - vy[dst]*dt)/l;
	
	int SX = sx;
	int SY = sy;
	
	TileId a00 = At<false>(SX, SY);
	TileId a10 = At<false>(SX+1, SY);
	TileId a11 = At<false>(SX+1, SY+1);
	TileId a01 = At<false>(SX, SY+1);
	
	// boundary check - to test
	{
		if(a00 < 0 && a10 < 0) {
			a00 = a01;
			a10 = a11;
		}
		if(a00 < 0 && a01 < 0) {
			a00 = a10;
			a01 = a11;
		}
		if(a11 < 0 && a10 < 0) {
			a11 = a01;
			a10 = a00;
		}
		if(a11 < 0 && a01 < 0) {
			a11 = a10;
			a01 = a00;
		}
		
		if(a00 < 0 || a01 < 0 || a10 < 0 || a11 < 0) {
			newSediment[dst] = 0;
			return;
		}
	}
	
	const float rx = sx - (float)SX;
	const float ry = sy - (float)SY;
	
	float ds00 = s[a00] * (1-rx) * (1-ry);
	float ds10 = s[a10] * (rx) * (1-ry);
	float ds01 = s[a01] * (1-rx) * (ry);
	float ds11 = s[a11] * (rx) * (ry);
	
	s[a00] -= ds00;
	s[a10] -= ds10;
	s[a01] -= ds01;
	s[a11] -= ds11;
	
	newSediment[dst] += ds00 + ds10 + ds01 + ds11;
}

template<bool safe>
void HydroErosion::SedimentTransportationUpdateWhatsLeft(int x, int y) {
	TileId src = At<false>(x, y);
	s[src] += newSediment[src];
	newSediment[src] = 0;
}




float HydroErosion::EvaporationRate(int x, int y) {
	return 0.001;
}

// 3.5
template<bool safe>
void HydroErosion::Evaporation(int x, int y) {
	TileId src = At<false>(x, y);
	d[src] *= 1 - EvaporationRate(x, y) * dt;
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
		return three{sumWater, sumSed, sumGround};
	};
	
	auto PrintWater = [&]() {
		auto w = SumWater();
// 		printf("water = %f  %f  %f   + %f\n", w.water, w.sediment, w.ground, w.sediment+w.ground);
	};
		
	
// 	printf("\n\n");
	PrintWater();
 	FOR_EACH_SAFE_BORDERS(1, CalcOutflux);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS(1, UpdateWaterLevelAndVelocity);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS(1, ErosionAndDeposition);
	PrintWater();
 	FOR_EACH_SAFE_BORDERS(0, SedimentTransportation);
 	FOR_EACH_SAFE_BORDERS(0, SedimentTransportationUpdateWhatsLeft);
	PrintWater();
	
// 	std::swap(newSediment, sediment);
	
 	FOR_EACH_SAFE_BORDERS(0, Evaporation);
}





