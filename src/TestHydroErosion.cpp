
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>

#include "lodepng.h"

#include "HydroErosion.hpp"

#include "File.hpp"

int main(int argc, char ** argv) {
	srand(time(NULL));
	
	float* data = NULL;
	unsigned width, height;
// 	Load(&data, width, height, 0, 100, argv[1]);
	width = height = 256;
	HydroErosion grid(width, height, 5);
// 	Grid grid;
	
	char str[1024];
	sprintf(str, "images/%s", argc > 1 ? argv[1] : "default");
	
	
// 	printf("Loaded [%u x %u]\n", width, height);
	data = new float[width*height];
// 	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
// 		grid.tiles[i].ground = data[i];
	for(int x=0; x<width; ++x) {
		for(int y=0; y<height; ++y) {
			float X = x/30.0f;
			float Y = y/30.0f;
			grid.ground[grid.At<false>(x, y)] = sin(X) * sin(Y)*100 * (
					sin(x/(float)(width-1) * 3.141592f) *
					sin(y/(float)(height-1) * 3.141592f)
					);
		}
	}
	uint32_t w, h;
	Load(grid.ground, w, h, 0, 1000, "../Noise256.png");
	printf("Converted\n");
	
	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
		data[i] = grid.ground[i];
	Save(data, width, height,
			(std::string(str)
			 + ".000.startpoint.png").c_str());
	
	grid.dt = 0.01;
	
	long long beg = clock();
	for(size_t I=0;; ++I) {
		for(int x=0; x<width; ++x) {
			for(int y=0; y<height; ++y) {
				grid.oldWater[grid.At<false>(x, y)] += 0.001 * grid.dt;
			}
		}
// 		for(int i=0; i<width*10; ++i) {
// 			grid.water[grid.At<false>(rand()%grid.width, rand()%grid.height)] += 0.01 * grid.dt;
// 		}
// 		grid.water[grid.At<false>(1330%width, 1850%height)] += 0.1 * grid.dt;
		grid.FullCycle();
		if(clock() - beg >= CLOCKS_PER_SEC*2) {
			printf(" done: %li ...", I);
			for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
				data[i] = grid.ground[i];
			Save(data, width, height,
					(std::string(str)
					 + "."
					 + std::to_string(I)
					 + ".eroding.png").c_str());
			printf(" saved\n");
			beg = clock();
		}
	}
	
	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
		data[i] = grid.ground[i];
	
	Save(data, width, height, (std::string(str) + ".eroded.png").c_str());
	delete[] data;
	return 0;
}

