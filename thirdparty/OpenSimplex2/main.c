#include <stdlib.h>
#include <stdio.h>
// #include "OpenSimplex2F/OpenSimplex2F.c"
#include "OpenSimplex2S/OpenSimplex2S.c"


#define BYTES_PER_PIXEL 3
#define FILE_HEADER_SIZE 14
#define INFO_HEADER_SIZE 40



unsigned char *createBitmapInfoHeader(int height, int width){
    static unsigned char infoHeader[] = {
        0,0,0,0, /// header size
        0,0,0,0, /// image width
        0,0,0,0, /// image height
        0,0,     /// number of color planes
        0,0,     /// bits per pixel
        0,0,0,0, /// compression
        0,0,0,0, /// image size
        0,0,0,0, /// horizontal resolution
        0,0,0,0, /// vertical resolution
        0,0,0,0, /// colors in color table
        0,0,0,0, /// important color count
    };

    infoHeader[ 0] = (unsigned char)(INFO_HEADER_SIZE);
    infoHeader[ 4] = (unsigned char)(width      );
    infoHeader[ 5] = (unsigned char)(width >>  8);
    infoHeader[ 6] = (unsigned char)(width >> 16);
    infoHeader[ 7] = (unsigned char)(width >> 24);
    infoHeader[ 8] = (unsigned char)(height      );
    infoHeader[ 9] = (unsigned char)(height >>  8);
    infoHeader[10] = (unsigned char)(height >> 16);
    infoHeader[11] = (unsigned char)(height >> 24);
    infoHeader[12] = (unsigned char)(1);
    infoHeader[14] = (unsigned char)(BYTES_PER_PIXEL*8);

    return infoHeader;
}

unsigned char *createBitmapFileHeader(int height, int stride){
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
        0,0,     /// signature
        0,0,0,0, /// image file size in bytes
        0,0,0,0, /// reserved
        0,0,0,0, /// start of pixel array
    };

    fileHeader[ 0] = (unsigned char)('B');
    fileHeader[ 1] = (unsigned char)('M');
    fileHeader[ 2] = (unsigned char)(fileSize      );
    fileHeader[ 3] = (unsigned char)(fileSize >>  8);
    fileHeader[ 4] = (unsigned char)(fileSize >> 16);
    fileHeader[ 5] = (unsigned char)(fileSize >> 24);
    fileHeader[10] = (unsigned char)(FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

void generateBitmapImage(unsigned char *image, int height, int width, char *imageFileName){
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;

    int stride = (widthInBytes) + paddingSize;

    FILE* imageFile = fopen(imageFileName, "wb");

    unsigned char* fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char* infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + (i*widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

void save_bitmap(char *filename, int width, int height, double **vals){
    unsigned char *image = (unsigned char *) malloc(sizeof(unsigned char) * height * width * BYTES_PER_PIXEL);
    int index = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            double val = vals[i][j];
//             val = (val + 1.0) / 2.0;
            unsigned char gray = (unsigned char) (val * 255);
            image[index++] = gray;
            image[index++] = gray;
            image[index++] = gray;
        }
    }

    generateBitmapImage(image, height, width, filename);
    free(image);
}



#define WIDTH 1024
#define HEIGHT 1024
#define PERIOD 128.0
#define OFF_X 2048
#define OFF_Y 2048	
#define FREQ 1.0 / PERIOD

OpenSimplexEnv *ose;
OpenSimplexGradients *osg;

OpenSimplexEnv *ose2;
OpenSimplexGradients *osg2;

OpenSimplexEnv *ose3;
OpenSimplexGradients *osg3;


float Noise(float x, float y, bool useNoise2)
{
	if (useNoise2) {
	return (noise2(ose, osg, x, y) + noise2(ose, osg, -x+100, -y+10) + 2.0) * 0.25;
	} else {
	return (noise2(ose, osg, x, y) + 1.0) * 0.5;
	}
}

float Height(float x, float y, float attenuation, int octaves, float coordMultAccu, bool useGrad, bool useNoise2, float heightScale)
{
	x = (x + OFF_X) * FREQ;
	y = (y + OFF_Y) * FREQ;
	
	float gx = 0, gy = 0;
	float dx = 0.001;
	
	float h = 0;
	float a = 1;
	float sum = 0;
	for (int i=0; i<octaves; ++i) {
		float h0, h1, h2;
		h0 = Noise(x, y, useNoise2) * a * heightScale;
		if (useGrad) {
			h1 = Noise(x+dx, y, useNoise2) * a * heightScale;
			h2 = Noise(x, y+dx, useNoise2) * a * heightScale;
			gx += (h1 - h0) / dx;
			gy += (h2 - h0) / dx;
			h += h0 / (1.0f + gx*gx + gy*gy);
		} else {
			h += h0;
		}
		
		sum += a;
		x *= coordMultAccu;
		y *= coordMultAccu;
		a *= attenuation;
	}
	return h / sum;
}

float AlternateBiomeNoise(float x, float y)
{
	x = (x + OFF_X) * FREQ;
	y = (y + OFF_Y) * FREQ;
	return noise2(ose3, osg3, x, y);
}

float ScaleNoise(float x, float y)
{
	x = (x + OFF_X) * FREQ;
	y = (y + OFF_Y) * FREQ;
	return noise2(ose2, osg2, -x-103, y) * 0.5 + 0.5;
}

float Terrain(float x, float y)
{
// 	return Height(x, y, 0.5, 1, 2.3, false, true);
	float biome = AlternateBiomeNoise(-x * 0.2, -y * 0.2);
	float scale = ScaleNoise(x*0.4, y*0.4);
	float mountains = Height(x, y, 0.5, 8, 2.3, true, true, scale * 0.6 + 0.4);
	float plains = Height(x+1000, -y+1000, 0.5, 8, 2.3, false, true, scale * 0.8 + 0.2) * 0.3;
	float h;
	if (biome < -0.3) {
		h = plains;
	} else if (biome < 0.3) {
		float f = (biome + 0.3) / 0.6;
		f = f * f * f * (f * (f * 6.0 - 15.0) + 10.0);
		h = plains * (1.0 - f) + f * mountains;
	} else{
		h = mountains;
	}
	return h;
}

int main(){
    ose = initOpenSimplex();
    osg = newOpenSimplexGradients(ose, 1234);
	
    ose2 = initOpenSimplex();
    osg2 = newOpenSimplexGradients(ose, 543249);
	
    ose3 = initOpenSimplex();
    osg3 = newOpenSimplexGradients(ose, 822199);
	
    double **noise = (double **) malloc(sizeof(double *) * HEIGHT);
    for (int y = 0; y < HEIGHT; y++){
        noise[y] = (double *) malloc(sizeof(double) * WIDTH);
        for (int x = 0; x < WIDTH; x++){
            noise[y][x] = Terrain(x, y-300);//Height(x, y, 0.5, 8, 2.3);
// 				(noise2(ose, osg, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ)
// 					+ noise2(ose2, osg2, (x + OFF_X) * FREQ, (y + OFF_Y) * FREQ)) * 0.5;
        }
    }
    save_bitmap("noise2.bmp", WIDTH, HEIGHT, noise);
    free(noise);
    return 0;
}
