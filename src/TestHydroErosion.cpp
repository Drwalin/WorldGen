
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>

#include <thread>
#include <chrono>

#include "lodepng.h"

#include "HydroErosion.hpp"

#include "File.hpp"

char str[1024];
void fileSaver(HydroErosion* grid, size_t* I) {
	size_t prevSaved = 0;
	while(true) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		if(*I == prevSaved)
			continue;
		printf(" done: %li ...", *I);
		prevSaved = *I;
		Save(grid->ground, grid->width, grid->height,
				(std::string(str)
				 + "."
				 + std::to_string(*I)
				 + ".eroding.png").c_str());
		printf(" saved\n");
	}
}

int main(int argc, char ** argv) {
	srand(time(NULL));
	
	float* data = NULL;
	unsigned width, height;
// 	Load(&data, width, height, 0, 100, argv[1]);
	width = height = 4096;
	HydroErosion grid(width, height);
// 	Grid grid;
	
	sprintf(str, "images/%s", argc > 1 ? argv[1] : "default");
	
	
// 	printf("Loaded [%u x %u]\n", width, height);
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
	Load(grid.ground, width, height, 0, 1000, "../PerlinSimplex.png");
	printf(" %f %f %f %f\n",
			grid.ground[grid.At<false>(0, height/2)],
			grid.ground[grid.At<false>(width-1, height/2)],
			grid.ground[grid.At<false>(width/2, 0)],
			grid.ground[grid.At<false>(width/2, height-1)]
		  );
	
	grid.dt = 0.1;
	
	size_t I=0;
	std::thread saverThread(fileSaver, &grid, &I);
	long long beg = clock();
	for(;; ++I) {
		for(int x=0; x<width; ++x) {
			for(int y=0; y<height; ++y) {
				grid.water[grid.At<false>(x, y)] += 0.01 * grid.dt;
			}
		}
// 		for(int i=0; i<width*10; ++i) {
// 			grid.water[grid.At<false>(rand()%grid.width, rand()%grid.height)] += 0.01 * grid.dt;
// 		}
// 		grid.water[grid.At<false>(1330%width, 1850%height)] += 0.1 * grid.dt;
		if(clock() - beg >= CLOCKS_PER_SEC) {
			grid.FullCycle(true);
			beg = clock();
		} else {
			grid.FullCycle(false);
		}
	}
	
	for(size_t i = 0; i<(size_t)width * (size_t)height; ++i)
		data[i] = grid.ground[i];
	
	Save(data, width, height, (std::string(str) + ".eroded.png").c_str());
	delete[] data;
	return 0;
}

