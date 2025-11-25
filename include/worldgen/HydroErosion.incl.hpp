#pragma once

#include <cassert>

#ifndef HYDRO_EROSION_INCL_HPP
#define HYDRO_EROSION_INCL_HPP

#include "HydroErosion.hpp"

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

#endif

