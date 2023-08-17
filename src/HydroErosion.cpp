
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
		int id = ((x+width)%width) + ((y+height)%height) * width;
		return id;
	} else {
		throw "HydroErosion::wrap == false is not supported as of now.";
		if constexpr(safe) {
			if(x < 0)            return -1;
			else if(x >= width)  return -1;
			if(y < 0)            return -1;
			else if(y >= height) return -1;
		}
		return x + y * width;
	}
}

template TileId HydroErosion::At<true>(int x, int y) const;
template TileId HydroErosion::At<false>(int x, int y) const;

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
	throw "Invalid direction given for TileId HydroErosion::Neighbour(int x, int y, int dir) const.";
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
float HydroErosion::GetMass(TileId id) const {
	return (oldWater[id] + oldSediment[id]) * A;
}

template<bool safe>
float HydroErosion::GetTotalHeight(TileId id) const {
	return ground[id] + oldWater[id] + oldSediment[id];
}

template<bool safe>
float HydroErosion::GetTotalHeight(glm::vec2 pos) const {
	throw "HydroErosion::GetTotalHeight(glm::vec2 pos) is unimplemented.";
	int x, y;
	x = pos.x;
	y = pos.y;
	float a00 = GetTotalHeight<safe>(At<safe>(x, y));
	float a01 = GetTotalHeight<safe>(At<safe>(x, y+1));
	float a10 = GetTotalHeight<safe>(At<safe>(x+1, y));
	float a11 = GetTotalHeight<safe>(At<safe>(x+1, y+1));
	
	const float rx = pos.x - (float)x;
	const float ry = pos.y - (float)y;
	
	float f00 = a00 * (1-rx) * (1-ry);
	float f01 = a01 * (1-rx) * (ry);
	float f10 = a10 * (rx) * (1-ry);
	float f11 = a11 * (rx) * (ry);
	
	return f00 + f01 + f10 + f11;
}





template<bool safe>
void HydroErosion::UpdateMomentun(int x, int y) {
	TileId src = At<safe>(x, y);
	TileId neighs[4];
	GetNeighboursOrSelf<safe>(x, y, src, neighs);
	
	const glm::vec2 dH = {
		GetTotalHeight<safe>(neighs[2])
			- GetTotalHeight<safe>(neighs[0]),
		GetTotalHeight<safe>(neighs[3])
			- GetTotalHeight<safe>(neighs[1]) };
	
	
	const glm::vec2 v = GetVelocity(src);
	const glm::vec2 delta = v*v + 2.f * g * (dH*0.5f);
	glm::vec2 t;
	if(delta.x < 0) {
		t.x = dt;
	} else {
		t.x = std::clamp<float>((-v.x + sqrt(delta.x)) / g, 0.f, dt);
	}
	if(delta.y < 0) {
		t.y = dt;
	} else {
		t.y = std::clamp<float>((-v.y + sqrt(delta.y)) / g, 0.f, dt);
	}
	
	glm::vec2 m = glm::min(dH*0.5f, oldWater[src] + oldSediment[src]);
	glm::vec2 dM = m * g * t;
	oldMomentum[src] += dM;
	
	
	
	
	
	
// 	glm::vec2 t = glm::clamp(glm::sqrt((glm::abs(dH) * 0.5f * 0.5f) * 2.f / g), 0.f, dt);
// 	glm::vec2 m = glm::min(dH*0.5f, oldWater[src] + oldSediment[src]);
// 	glm::vec2 dM = m * g * t;
// 	oldMomentum[src] += dM;
// 	printf("dM = %f, %f\n", dM.x, dM.y);
	
	
	
	
	float len = glm::length(oldMomentum[src]);
	
// 	oldMomentum[src] += (dH * 0.5f * g * dt * (oldWater[src] + oldSediment[src])) / (len/10.f+1.f);
	oldMomentum[src] *= 0.98f;
	
	if(len > 100) {
		oldMomentum[src] *= 100.f/len;
	}
}

glm::vec2 HydroErosion::GetVelocity(TileId id) const {
	float mass = GetMass<false>(id);
	if(mass < 0.001) {
		oldMomentum[id] = {0,0};
		return {0,0};
	}
	return oldMomentum[id]/mass;
}



template<bool safe>
float HydroErosion::GetSedimentCapacity(TileId t, int x, int y) const {
	TileId neighs[4];
	GetNeighboursOrSelf<safe>(x, y, neighs);
	float h[4];
	for(int i=0; i<4; ++i) {
		h[i] = ground[neighs[i]];
	}
	float dbdx = (h[2]-h[0]) / (2.f * l);
	float dbdy = (h[3]-h[1]) / (2.f * l);
	float dxy = dbdx*dbdx + dbdy*dbdy;
	float sin_tilt = sqrt(dxy) / sqrt(1 + dxy);
	return std::max(1.f,
			Kc * (sin_tilt + 0.15f)
				* std::max(0.15f,
						glm::length(GetVelocity(t)))
				* oldWater[t]);
}

template<bool safe>
void HydroErosion::UpdateErosionAndDeposition(int x, int y) {
	TileId src = At<safe>(x, y);
	float capacity = GetSedimentCapacity<safe>(src, x, y);
	if(capacity < 0.001) {
		capacity = 0.001;
	}
	if(capacity < oldSediment[src]) {
		const float dm = (oldSediment[src]-capacity) * Ks * dt;
		const float prevMass = oldSediment[src] + oldWater[src];
		const float momentumFactor = (prevMass - dm) / prevMass;
		oldMomentum[src] *= momentumFactor;
		ground[src] += dm;
		oldSediment[src] -= dm;
	} else {
		const float dm = (capacity-oldSediment[src]) * Kd * dt;
		ground[src] -= dm;
		oldSediment[src] += dm;
	}
}





template<bool safe>
void HydroErosion::SedimentWaterAndMomentumTransportation(int x, int y) {
	TileId src = At<false>(x, y);
	const glm::vec2 v = GetVelocity(src);
	const float sx = (x*l + v.x*dt)/l;
	const float sy = (y*l + v.y*dt)/l;
	
	int SX = floor(sx);
	int SY = floor(sy);
	
	const TileId ids[4] = {
		At<false>(SX, SY),
		At<false>(SX+1, SY),
		At<false>(SX+1, SY+1),
		At<false>(SX, SY+1) };
	
	const float rx = sx - (float)SX;
	const float ry = sy - (float)SY;
	
	const float f[4] = {
		(1-rx) * (1-ry),
		(rx) * (1-ry),
		(rx) * (ry),
		(1-rx) * (ry) };
	
	float sediment = oldSediment[src];
	float water = oldWater[src];
	glm::vec2 momentum = oldMomentum[src];
	
	for(int i=0; i<4; ++i) {
		newSediment[ids[i]] += sediment * f[i];
		newWater[ids[i]] += water * f[i];
		newMomentum[ids[i]] += momentum * f[i];
	}
	
	oldWater[src] = 0;
	oldSediment[src] = 0;
	oldMomentum[src] = {0,0};
}





template<bool safe>
void HydroErosion::Evaporation(int x, int y) {
	TileId src = At<false>(x, y);
	float dw = std::max(std::min<float>(0.1, newWater[src] * 0.1), newWater[src]);
	
	float divider = (newWater[src] - dw + newSediment[src]);
	if(divider < 0.0001) {
		newMomentum[src] = {0,0};
	} else {
		float momentumFactor = (newWater[src] - dw + newSediment[src])
			/ divider;
		
		newMomentum[src] *= momentumFactor;
	}
	newWater[src] -= dw;
}





template<bool safe>
void HydroErosion::ThermalErosionStep1(int x, int y) {
	// angle 45*
	const float thermalErosionRate = 1;
	TileId src = At<safe>(x, y);
	
	TileId neighs[4];
	GetNeighboursOrSelf<safe>(x, y, neighs);
	float h[4];
	for(int i=0; i<4; ++i) {
		h[i] = ground[neighs[i]];
	}
	glm::vec2 dd;
	dd.x = (h[2]-h[0]) / (2.f * l);
	dd.y = (h[3]-h[1]) / (2.f * l);
	
	float l = glm::length(dd);
	if(l < 0.01) {
		int r = rand()%4;
		TileId neigh = Neighbour<safe>(x, y, r);
		if(ground[neigh]+l < ground[src]) {
			float dh = ((ground[src] - (ground[neigh]+l))*0.5f) * thermalErosionRate * dt;
			ground[src] -= dh;
			ground[neigh] += dh;
		}
		return;
	}
	
	dd *= 0.2f/l;
	
	
	float dstH = 0;
	glm::vec2 pos = dd + glm::vec2{x,y};
	float rx, ry;
	{
		x = pos.x;
		y = pos.y;
		float a00 = ground[At<safe>(x, y)];
		float a01 = ground[At<safe>(x, y+1)];
		float a10 = ground[At<safe>(x+1, y)];
		float a11 = ground[At<safe>(x+1, y+1)];
		
		rx = pos.x - (float)x;
		ry = pos.y - (float)y;
		
		float f00 = a00 * (1-rx) * (1-ry);
		float f01 = a01 * (1-rx) * (ry);
		float f10 = a10 * (rx) * (1-ry);
		float f11 = a11 * (rx) * (ry);
		
		dstH = f00 + f01 + f10 + f11;
	}
	dstH += l;
	
	float srcH = ground[src];
	
	if(srcH > dstH + 1) {
		float dh = ((srcH - dstH)*0.5f) * thermalErosionRate * dt;
		ground[src] -= dh;
		
		ground[At<safe>(x, y)] += dh * (1-rx) * (1-ry);
		ground[At<safe>(x, y+1)] += dh * (1-rx) * (ry);
		ground[At<safe>(x+1, y)] += dh * (rx) * (1-ry);
		ground[At<safe>(x+1, y+1)] += dh * (rx) * (ry);
	}
}

template<bool safe>
void HydroErosion::ThermalErosionStep2(int x, int y) {
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
			sumWater += oldWater[i];
		}
		float sumGround = 0;
		for(int i=0; i<width*height; ++i) {
			sumGround += ground[i];
		}
		float sumSed = 0;
		for(int i=0; i<width*height; ++i) {
			sumSed += oldSediment[i];
		}
		return three{sumWater/(width*height), sumSed/(width*height), sumGround/(width*height)};
	};
	
	auto PrintWater = [&]() {
		auto w = SumWater();
		printf("water: %f,  sediment: %f,  ground: %f   all dirt: %f\n", w.water, w.sediment, w.ground, w.sediment+w.ground);
	};
		
	
// 	printf("\n\n");
// 	PrintWater();
 	FOR_EACH_SAFE_BORDERS(0, UpdateMomentun);
// 	PrintWater();
 	FOR_EACH_SAFE_BORDERS(0, UpdateErosionAndDeposition);
// 	PrintWater();
 	FOR_EACH_SAFE_BORDERS(0, SedimentWaterAndMomentumTransportation);
	
 	FOR_EACH_SAFE_BORDERS(0, Evaporation);
	
	std::swap(oldWater, newWater);
	std::swap(oldMomentum, newMomentum);
	std::swap(oldSediment, newSediment);
	
 	FOR_EACH_SAFE_BORDERS(0, ThermalErosionStep1);
 	FOR_EACH_SAFE_BORDERS(0, ThermalErosionStep2);
}





/*





template<bool safe, int dir>
float HydroErosion::CalcFluxInDirection(TileId src, TileId neigh) const {
	const float f_old = fluxArray[src].v[dir];
	const float dh = b[src] + d[src] - b[neigh] - d[neigh];
// 	const float dh = b[src] + s[src] + d[src] - b[neigh] - s[neigh] - d[neigh];
	return std::max<float>(
			0.f,
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
	fluxArray[src].v[0] = 0;
	fluxArray[src].v[1] = 0;
	fluxArray[src].v[2] = 0;
	fluxArray[src].v[3] = 0;
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
	const float sin_tilt = std::max<float>(SinusLocalTiltAngle<safe>(src, x, y), 0.0);
	float C = Kc * sin_tilt * sqrt(vx[src]*vx[src] + vy[src]*vy[src]);
// 	C = std::min(C, 2.0f);
// 	printf(" C = %f (v = %f, %f)\n", C, vx[src], vy[src]);
	
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
// 	TileId dst = At<false>(x, y);
// 	float sx, sy;
// 	sx = (x - vx[dst]*dt)/l;
// 	sy = (y - vy[dst]*dt)/l;
// 	
// 	int SX = sx;
// 	int SY = sy;
// 	
// 	TileId a00 = At<false>(SX, SY);
// 	TileId a10 = At<false>(SX+1, SY);
// 	TileId a11 = At<false>(SX+1, SY+1);
// 	TileId a01 = At<false>(SX, SY+1);
// 	
// 	// boundary check - to test
// 	{
// 		if(a00 < 0 && a10 < 0) {
// 			a00 = a01;
// 			a10 = a11;
// 		}
// 		if(a00 < 0 && a01 < 0) {
// 			a00 = a10;
// 			a01 = a11;
// 		}
// 		if(a11 < 0 && a10 < 0) {
// 			a11 = a01;
// 			a10 = a00;
// 		}
// 		if(a11 < 0 && a01 < 0) {
// 			a11 = a10;
// 			a01 = a00;
// 		}
// 		
// 		if(a00 < 0 || a01 < 0 || a10 < 0 || a11 < 0) {
// 			newSediment[dst] = 0;
// 			return;
// 		}
// 	}
// 	
// 	const float rx = sx - (float)SX;
// 	const float ry = sy - (float)SY;
// 	
// 	float ds00 = s[a00] * (1-rx) * (1-ry);
// 	float ds10 = s[a10] * (rx) * (1-ry);
// 	float ds01 = s[a01] * (1-rx) * (ry);
// 	float ds11 = s[a11] * (rx) * (ry);
// 	
// 	s[a00] -= ds00;
// 	s[a10] -= ds10;
// 	s[a01] -= ds01;
// 	s[a11] -= ds11;
// 	
// 	newSediment[dst] += ds00 + ds10 + ds01 + ds11;
	
	
	
	TileId dst = At<false>(x, y);
	float sx, sy;
	sx = (x + vx[dst]*dt)/l;
	sy = (y + vy[dst]*dt)/l;
	
	int SX = sx;
	int SY = sy;
	
	TileId a00 = At<true>(SX, SY);
	TileId a10 = At<true>(SX+1, SY);
	TileId a11 = At<true>(SX+1, SY+1);
	TileId a01 = At<true>(SX, SY+1);
	
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
			newSediment[dst] = s[dst];
			return;
		}
	}
	
	const float rx = sx - (float)SX;
	const float ry = sy - (float)SY;
	
	const float sed = s[dst];
	s[dst] = 0;
	
	float ds00 = sed * (1-rx) * (1-ry);
	float ds10 = sed * (rx) * (1-ry);
	float ds01 = sed * (1-rx) * (ry);
	float ds11 = sed * (rx) * (ry);
	
	newSediment[a00] += ds00;
	newSediment[a10] += ds10;
	newSediment[a01] += ds01;
	newSediment[a11] += ds11;
}

template<bool safe>
void HydroErosion::SedimentTransportationUpdateWhatsLeft(int x, int y) {
	TileId src = At<false>(x, y);
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
	d[src] -= std::min<float>(0.1, d[src] * EvaporationRate(x, y)) * dt;
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
		
	
// 	printf("\n\n");
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


*/


