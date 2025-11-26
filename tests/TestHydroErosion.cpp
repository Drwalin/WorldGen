
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "../thirdparty/lodepng/lodepng.h"

#include "../include/worldgen/HydroErosion.hpp"

#include "../include/worldgen/File.hpp"

int main(int argc, char ** argv) {
	srand(time(NULL));
	Grid grid;
	
	char str[1024];
	sprintf(str, "images/%s", argc > 1 ? argv[1] : "default");
	
	float* data = NULL;
	unsigned width, height;
// 	Load(&data, width, height, 0, 100, argv[1]);
	width = height = 256;
	
	
	printf("Loaded [%u x %u]\n", width, height);
	grid.Init(width, height);
	data = new float[width*height];
// 	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
// 		grid.tiles[i].ground = data[i];
	for(int x=0; x<width; ++x) {
		for(int y=0; y<height; ++y) {
			float X = x/30.0f;
			float Y = y/30.0f;
			grid.ground[grid.At<false>(x, y)] = sin(X) * sin(Y)*100;
		}
	}
	printf("Converted\n");
	
	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
		data[i] = grid.ground[i];
	Save(data, width, height, "Generated.png");
	
	long long beg = clock();
	for(size_t I=0;; ++I) {
		for(int x=0; x<width; ++x) {
			for(int y=0; y<height; ++y) {
				grid.water[grid.At<false>(x, y)] += 0.01;
			}
		}
// 		for(int i=0; i<width*10; ++i) {
// 			grid.At<false>(rand()%grid.width, rand()%grid.height)->water += 0.01;
// 		}
// 		grid.At<false>(1330%width, 1850%height)->water += 0.1;
		grid.FullCycle();
// 		if(clock() - beg >= CLOCKS_PER_SEC*2) {
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
// 		}
	}
	
	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
		data[i] = grid.ground[i];
	
	Save(data, width, height, (std::string(str) + ".eroded.png").c_str());
	delete[] data;
	return 0;
}

