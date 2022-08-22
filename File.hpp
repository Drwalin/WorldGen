
#include <vector>
#include <cinttypes>

#include "lodepng.h"

void Save(float* data, unsigned width, unsigned height, const char* filename) {
	float min = data[0], max = data[0];
	for(size_t i=0; i<(size_t)width*(size_t)height; ++i) {
		if(min > data[i])
			min = data[i];
		else if(max < data[i])
			max = data[i];
	}
	uint16_t* image = new uint16_t[width*height];
	float m = ((float)(0xFFFF))/(max-min);
	for(size_t i=0; i<(size_t)width*(size_t)height; ++i) {
		int v = (data[i] - min) * m;
		if(v < 0)
			v = 0;
		else if(v > 0xFFFF)
			v = 0xFFFF;
		image[i] = v;
	}
	printf(" min:max = %f:%f\n", min, max);
	lodepng_encode_file(filename, (const unsigned char*)image, width ,height,
			LCT_GREY, 16);
	delete[] image;
}

void Load(float** data, unsigned& width, unsigned &height, float min, float max,
		const char* filename) {
	unsigned char* image = NULL;
	lodepng_decode_file(&image, &width, &width, filename, LCT_GREY, 16);
	uint16_t* ptr = (uint16_t*)image;
	*data = new float[width*height];
	for(size_t i=0; i<(size_t)width*(size_t)height; ++i) {
		(*data)[i] = (((float)(ptr[i])) / ((float)0xFFFF)) * (max - min) + min;
	}
}

