#include <chrono>
#include <thread>
#include <atomic>
#include <random>

#include <glm/fwd.hpp>
#include <glm/geometric.hpp>
#include <glm/vector_relational.hpp>

#include "../../../OpenGLWrapper/samples/Camera.hpp"
#include "../../../OpenGLWrapper/samples/DefaultCameraAndOtherConfig.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/VBO.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/VAO.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/OpenGL.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Shader.hpp"
#include "../../../OpenGLWrapper/include/openglwrapper/Texture.hpp"
// #include "../../../OpenGLWrapper/include/openglwrapper/Sync.hpp"
// #include "../../../OpenGLWrapper/include/openglwrapper/BufferAccessor.hpp"

#include "../../include/worldgen/HydroErosion.hpp"
#include "../../include/worldgen/Noises.hpp"

volatile bool useGpu = true;
volatile bool gridInited = false;
int width = 512 + 64;
int height = width;
int riverSourcesCount = 4;

float generatorYScale = 1.0;
float SOFT_LAYER = 0.1;

float uniScale = 0.01;

float SCALE = 10;

float HYDRO_EROSION_Y_SCALE = 1.0f;

float horizontalScale = 1.0f * uniScale;
float verticalScale = uniScale;//1200.0f * uniScale / SCALE;
float noiseHorizontalScale = 1.0f / 4000.0f * SCALE;
volatile int HYDRO_ITER = 0;
volatile float averageHydroIterationDuration = 0.0f;
volatile bool disableSimulation = false;
volatile bool updateColorsTexture = false;

volatile long double SUM_MATERIAL = 0;

std::mt19937_64 mt(732168281);

struct VertexHeights {
	float h = 0;
	float w = 0;
};
VertexHeights *vertHeights = nullptr;

struct VertexColors {
	uint8_t rgba[4] = {0, 0, 0, 1};
};
VertexColors *vertColors = nullptr;

void ThreadFunction();
void HydroErosionIteration();

volatile float hydroErosionDuration = 0;
volatile bool updateHeights = true;
volatile bool updateWaterHeights = true;

float clamp(float v, float min, float max)
{
	if (v < min) {
		return v;
	} else if (v > max) {
		return max;
	}
	return v;
}

int clamp(int v, int min, int max)
{
	if (v < min) {
		return v;
	} else if (v > max) {
		return max;
	}
	return v;
}

Grid grid;

int main(int argc, char **argv)
{
	printf("sizeof(glm::vec3) = %lu\n", sizeof(glm::vec3));
	
	width = argc > 2 ? atoi(argv[2]) : 512 + 64;
	if (width < 512 + 64) {
		width = 512+64;
	}
	if (width > 16384) {
		width = 16384;
	}
	height = width;
	riverSourcesCount = argc > 3 ? clamp(atoi(argv[3]), 0, 1000) : 4;
	generatorYScale = argc > 4 ? clamp((float)atof(argv[4]), 0.001f, 100.0f) : 1.0f;
	SOFT_LAYER = argc > 5 ? clamp((float)atof(argv[4]), 0.0f, 1000.0f) : 1.0f;
	const int maxMeshSize = argc > 1 ? atoi(argv[1]) : 512;
	const int meshWidth = width > maxMeshSize ? maxMeshSize : width;
	const int meshHeight = height > maxMeshSize ? maxMeshSize : height;
	
	DefaultsSetup(false);

	// Load shader
	gl::Shader shader;
	shader.Load("../tests/gen_3d/vertex.glsl", "../tests/gen_3d/geometry.glsl",
				"../tests/gen_3d/fragment.glsl");

	// TODO: load from argv

	// Generate index data
	gl::VBO ebo(sizeof(uint32_t), gl::ELEMENT_ARRAY_BUFFER, gl::DYNAMIC_DRAW);
	ebo.Init();
	
	vertHeights = new VertexHeights[width*height]({0.0f, 0.0f});
	vertColors = new VertexColors[width*height]({0, 0, 0, 255});
	
	gl::Texture heightsTexture, colorsTexture;
	
// 	heightsTexture.WrapX(gl::TextureWrapParam::CLAMP_TOtEDGE);
// 	heightsTexture.WrapY(gl::TextureWrapParam::CLAMP_TOtEDGE);
// 	heightsTexture.MagFilter(gl::TextureMagFilter::MAG_LINEAR);
// 	heightsTexture.MinFilter(gl::TextureMinFilter::LINEAR);
// 	
	heightsTexture.InitTextureEmpty(width, height, gl::TextureTarget::TEXTURE_2D, gl::TextureSizedInternalFormat::RG32F);
	colorsTexture.InitTextureEmpty(width, height, gl::TextureTarget::TEXTURE_2D, gl::TextureSizedInternalFormat::RGBA8);
	
// 	heightsTexture.WrapX(gl::TextureWrapParam::CLAMP_TOtEDGE);
// 	heightsTexture.WrapY(gl::TextureWrapParam::CLAMP_TOtEDGE);
// 	heightsTexture.MagFilter(gl::TextureMagFilter::MAG_LINEAR);
// 	heightsTexture.MinFilter(gl::TextureMinFilter::LINEAR);
// 	
	colorsTexture.Update2((const void*)vertColors, 0, 0, width, height, 0, gl::TextureDataFormat::RGBA, gl::DataType::UNSIGNED_BYTE);
	heightsTexture.Update2((const void*)vertHeights, 0, 0, width, height, 0, gl::TextureDataFormat::RG, gl::DataType::FLOAT);
	
	bool filterWrapTextureSet = false;
	
	heightsTexture.WrapX(gl::TextureWrapParam::CLAMP_TOtEDGE);
	heightsTexture.WrapY(gl::TextureWrapParam::CLAMP_TOtEDGE);
	heightsTexture.MagFilter(gl::TextureMagFilter::MAG_LINEAR);
	heightsTexture.MinFilter(gl::TextureMinFilter::LINEAR);

	{
		std::vector<uint32_t> Ebo;
		Ebo.reserve(3 * 2 * (meshWidth - 1) * (meshHeight - 1));
		int yStride = 8;
		for (int Y = 0; Y + 1 < meshHeight; Y+=yStride) {
			for (int x = 0; x + 1 < meshWidth; ++x) {
				for (int y = Y; y + 1 < meshHeight && y < Y+yStride; ++y) {
					int id0 = x + y * meshWidth;
					int id1 = id0 + 1;
					int id2 = id0 + meshWidth;
					int id3 = id1 + meshWidth;
					Ebo.push_back(id0);
					Ebo.push_back(id1);
					Ebo.push_back(id2);
					Ebo.push_back(id2);
					Ebo.push_back(id1);
					Ebo.push_back(id3);
				}
			}
		}
		// Generate VBO from vertex data
		ebo.Generate(Ebo.data(), Ebo.size());
	}

	std::thread(ThreadFunction).detach();

	// Initiate VAO with VBO attributes
	gl::VAO vao(gl::TRIANGLES);
	vao.Init();
	vao.BindElementBuffer(ebo, gl::DataType::UNSIGNED_INT);

	// Get uniform locations
	const int modelLoc = shader.GetUniformLocation("model");
	const int viewLoc = shader.GetUniformLocation("view");
	const int projLoc = shader.GetUniformLocation("projection");
	const int useWaterLoc = shader.GetUniformLocation("useWater");
	const int gridWidthLoc = shader.GetUniformLocation("gridWidth");
	const int colorTexLoc = shader.GetUniformLocation("colorTex");
	const int heightTexLoc = shader.GetUniformLocation("heightTex");
	const int meshSizeLoc = shader.GetUniformLocation("meshSize");
	const int cameraPosLoc = shader.GetUniformLocation("cameraPos");
	const int centerOffsetLoc = shader.GetUniformLocation("centerOffset");
	const int generatorYScaleLoc = shader.GetUniformLocation("generatorYScale");
	
	shader.SetFloat(generatorYScaleLoc, generatorYScale);
	
	shader.SetTexture(colorTexLoc, &colorsTexture, 0);
	shader.SetTexture(heightTexLoc, &heightsTexture, 1);

	shader.Use();
	glUniform2iv(shader.GetUniformLocation("size"), 1,
				 std::vector<int32_t>{width, height}.data());
	glUniform2iv(meshSizeLoc, 1, std::vector<int32_t>{meshWidth, meshHeight}.data());
	GL_CHECK_PUSH_ERROR;
	shader.SetVec3(shader.GetUniformLocation("scale"),
				   glm::vec3(horizontalScale, verticalScale, horizontalScale));
	shader.SetInt(gridWidthLoc, width);

	// Init camera position
	camera.ProcessMouseMovement(175, 50);
	camera.position = glm::vec3(-25, 65, -25);

	auto beg = std::chrono::steady_clock::now();
	int frames = 0;
	float fps = 0;
	bool disableRender = false;
	
	auto lastUpdateHeightsTexture = std::chrono::steady_clock::now() - std::chrono::seconds(10);
	
	while (!glfwWindowShouldClose(gl::openGL.window)) {
		
		++frames;
		if (frames % 50 == 0) {
			auto now = std::chrono::steady_clock::now();
			auto dt = now - beg;
			auto ns = dt.count();
			auto sec = ns / 1'000'000'000.0;
			fps = frames / sec;
			beg = now;
			frames = 0;
		}
		printf("\r    fps: %5.2f         (%6i) hydro: %7.2f ms  (%7.2f ms)     mat = %Lf                            ", fps, HYDRO_ITER,
			   hydroErosionDuration, averageHydroIterationDuration, SUM_MATERIAL);
		fflush(stdout);

		DefaultIterationStart();
		
		if (gl::openGL.IsKeyDown('U')) {
			updateHeights = true;
		}
		if (gl::openGL.WasKeyPressed('I')) {
			updateWaterHeights ^= true;
		}
		
		if (gl::openGL.WasKeyPressed('O')) {
			disableRender = !disableRender;
		}
		if (gl::openGL.WasKeyPressed('P')) {
			disableSimulation = !disableSimulation;
		}

		{
			if (updateColorsTexture) {
				updateColorsTexture = false;
				colorsTexture.Update2((const void*)vertColors, 0, 0, width, height, 0, gl::TextureDataFormat::RGBA, gl::DataType::UNSIGNED_BYTE);
				
				gridInited = true;
			}
			
			if (gridInited && useGpu) {
				HydroErosionIteration();
			}
			
			if (updateWaterHeights || updateHeights) {
				const auto now = std::chrono::steady_clock::now();
				auto dur = (now - lastUpdateHeightsTexture);
				if (dur > std::chrono::milliseconds(500)) {
					lastUpdateHeightsTexture = now;
					if (gridInited && useGpu) {
						grid.gpu.UpdateHeightsTexture(&heightsTexture);
					} else {
						heightsTexture.Update2((const void*)vertHeights, 0, 0, width, height, 0, gl::TextureDataFormat::RG, gl::DataType::FLOAT);
					}
					if (filterWrapTextureSet == false) {
						filterWrapTextureSet = true;

						heightsTexture.WrapX(gl::TextureWrapParam::CLAMP_TOtEDGE);
						heightsTexture.WrapY(gl::TextureWrapParam::CLAMP_TOtEDGE);
						heightsTexture.MagFilter(gl::TextureMagFilter::MAG_LINEAR);
						heightsTexture.MinFilter(gl::TextureMinFilter::LINEAR);
					}
				}
			}
		}

		// Use shader
		shader.Use();

		// Calculate projection matrix
		glm::mat4 projection = glm::perspective(
			45.0f, (float)gl::openGL.GetWidth() / (float)gl::openGL.GetHeight(),
			0.1f, 10000.0f);

		// Calculate view matrix
		glm::mat4 view = camera.getViewMatrix();

		// Calulate model matrix
		glm::mat4 model =
			glm::translate(glm::scale(glm::mat4(1.0f), glm::vec3(10)),
						   glm::vec3(0.0f, -1.0f, 0.0f));

		// Set shader uniform matrices
		shader.SetMat4(viewLoc, view);
		shader.SetMat4(projLoc, projection);
		shader.SetMat4(modelLoc, model);
		{
			glm::vec3 cp = camera.position * 10.0f;
			shader.SetVec3(cameraPosLoc, camera.position * 10.0f);
			glm::ivec2 of = {cp.x, cp.z};
			if (of.x < meshWidth/2) {
				of.x = meshWidth/2;
			} else if (of.x > width-meshWidth/2 - 1) {
				of.x = width-meshWidth/2 - 1;
			}
			if (of.y < meshHeight/2) {
				of.y = meshHeight/2;
			} else if (of.y > height-meshHeight/2 - 1) {
				of.y = height-meshHeight/2 - 1;
			}
			shader.SetVec2(centerOffsetLoc, glm::vec2(of) / (glm::vec2(width, height)));
		}

		// Draw VAO
		if (disableRender == false) {
			shader.SetTexture(colorTexLoc, &colorsTexture, 0);
			shader.SetTexture(heightTexLoc, &heightsTexture, 1);
			shader.SetInt(useWaterLoc, 0);
			glDepthMask(true);
			vao.Draw();
			shader.SetInt(useWaterLoc, 1);
			glDepthMask(false);
			vao.Draw();
			glDepthMask(true);
		}

		DefaultIterationEnd();
		
		gl::openGL.PrintErrors();
		
		if (useGpu == false || gridInited == false || disableRender) {
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
	}

	gl::openGL.Destroy();
	glfwTerminate();

	return 0;
}

glm::ivec3 ColorGradient(int x, int y)
{
	thread_local wg::SimplexNoise colorSimplex(12342312);

	float n = colorSimplex.Noise2(glm::vec3{x, 0, y} * 130.0f);

	glm::vec3 c = {n, n, n};
	c = glm::vec3(n * 230.0f);
	c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
	return c;
}

void ThreadFunction()
{
	wg::SimplexNoise simplex(1234);

	int CHUNK_SIZE = 128;

	std::vector<glm::ivec2> chunks;

	for (int _x = 0; _x < width; _x += CHUNK_SIZE) {
		for (int _y = 0; _y < height; _y += CHUNK_SIZE) {
			chunks.push_back({_x, _y});
		}
	}

	std::vector<std::thread> threads;

	std::atomic<int> chunkCounter = 0;
	std::atomic<int> threadsDone = 0;

	auto threadCalcFunc = [&]() {
		float maxH = 0;
		for (;;) {
			int id = chunkCounter.fetch_add(1);
			if (id >= chunks.size()) {
				break;
			}
			glm::ivec2 chunk = chunks[id];
			for (int _y = chunk.y; _y < height && _y < chunk.y + CHUNK_SIZE;
				 ++_y) {
				for (int _x = chunk.x; _x < width && _x < chunk.x + CHUNK_SIZE;
					 ++_x) {
					glm::vec2 p{_x, _y};
// 					p -= glm::vec2(2000, 1000);
					p *= 0.6f;
					p += 2700.0f;
					p *= 0.5f;
					
					
					glm::vec3 v;

					p.y += 750.f / SCALE;

					p.x += (simplex.Fbm(glm::vec2(-p.x / 53 + 100, p.y / 53 - 1000),
									  3, 0.5, 2.3, false, false, 1.0f) -
						  0.5) *
						 10;
					p.y += (simplex.Fbm(glm::vec2(+p.x / 53 - 1000, -p.y / 53 + 100),
									  3, 0.5, 2.3, false, false, 1.0f) -
						  0.5) *
						 10;

					v.x =
						simplex.Terrain(glm::vec2{p.x, p.y} * noiseHorizontalScale,
										horizontalScale);

					maxH = std::max(maxH, v.x);

					float h = v.x * 600.0f * generatorYScale;
					glm::ivec3 c = ColorGradient(_x, _y);

					const int i = _x * height + _y;
					vertColors[i] = {(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 255};
					vertHeights[i] = {h, 0};
				}
			}
		}
		printf("Max height: %.3f\n", maxH);
		threadsDone++;
	};

	int threadsCount =
		std::max(((int)std::thread::hardware_concurrency()) - 2, 1);

	for (int i = 1; i < threadsCount; ++i) {
		threads.push_back(std::thread(threadCalcFunc));
		threads.back().detach();
	}

	threadCalcFunc();

	while (threadsDone.load() < threads.size() + 1) {
		std::this_thread::sleep_for(std::chrono::milliseconds(60));
	}
	printf("\r Done!                         \n");
	updateColorsTexture = true;

	if (useGpu == false) {
		for (;;) {
			HydroErosionIteration();
		}
	}
}

void HydroErosionIteration()
{
	static auto lastUpdateHeightsTexture = std::chrono::steady_clock::now() - std::chrono::seconds(10);
	static wg::SimplexNoise simplex(432127);
	static std::vector<glm::vec3> riverSources;
	
	if (HYDRO_ITER == 0) {
		grid.Init(width, height, useGpu);
		
		for (int i=0; i<riverSourcesCount; ++i) {
			glm::vec3 p;
			p.x = ((mt() % width) - 128) + 64;
			p.y = ((mt() % height) - 128) + 64;
			p.z = (mt() % 1000) / 100.0f + 10.0f;
			riverSources.push_back(p);
		}

		int maxHeight = 0;
		wg::SimplexNoise simplex(13222);
		for (int _y = 0; _y < height; ++_y) {
			for (int _x = 0; _x < width; ++_x) {
				const int i = _x + _y * width;
				int t = grid.At(_x, _y);
				grid.ground[t].layers[0] = vertHeights[i].h * HYDRO_EROSION_Y_SCALE - SOFT_LAYER;
				grid.ground[t].layers[1] = SOFT_LAYER;
				if (_x > 64 && _x < width - 64 && _y > 64 && _y < height - 64) {
					if (vertHeights[i].h > vertHeights[maxHeight].h) {
						maxHeight = i;
						if (riverSources.size() > 0) {
							riverSources[0] = {_x, _y, 500};
						}
					}
				}
				grid.water[t] = 0.1;
			}
		}
		if (riverSources.size() > 0) {
			riverSources[0].x += (mt() % 11) - 5;
			riverSources[0].y += (mt() % 11) - 5;
		}
		
		HYDRO_ITER = HYDRO_ITER + 1;
		
		if (useGpu) {
			grid.gpu.UpdateGround(grid.ground);
			grid.gpu.UpdateWater(grid.water);
			grid.gpu.UpdateSediment(grid.sediment);
			grid.gpu.UpdateTemp1(grid.temp1);
			grid.gpu.UpdateVelocity(grid.velocity);
			grid.gpu.UpdateFlux(grid.flux);
			grid.gpu.UpdateRiverSources(riverSources.data(), riverSources.size());
		}
	}
	
	if (disableSimulation) {
		if (useGpu == false) {
			std::this_thread::sleep_for(std::chrono::milliseconds(150));
		}
		return;
	}
	
	const auto a = std::chrono::steady_clock::now();
	
	if (HYDRO_ITER % 1000 == 0) {
		for (int _y = 0; _y < height; ++_y) {
			for (int _x = 0; _x < width; ++_x) {
				const float x = _x;
				const float y = _y;

				const float noise =
					simplex.Noise2(glm::vec3(-x / 531 - 100, y / 531 + 1000, HYDRO_ITER / 1000.0f));

				int p = grid.At(_x, _y);

				grid.water[p] += noise * 0.1;
// 					grid.ground[p].layers[0] += noise * 3.0f;
			}
		}
	}
	
	if (grid.useWater) {
		for (auto p : riverSources) {
			int i = grid.At(p.x, p.y);
			grid.water[i] += p.z * grid.dt;
		}
		
		for (int i=0; i<width * pow(width, 0.2); ++i) {
			int id = (mt()%(width*height));
			grid.water[id] += 0.01f;
		}

		if (HYDRO_ITER % 500 == 0) {
			long double SUM = 0;
			for (int _y = 0; _y < height; ++_y) {
				for (int _x = 0; _x < width; ++_x) {
					/*
					const float x = _x;
					const float y = _y;

					const float _rain =
						simplex.Noise2(glm::vec3(-x / 531 - 100, y / 531 + 1000, HYDRO_ITER / 100.0f))
						*2.0f - 0.5f;
					
					const float rain = _rain < 0.0f ? 0.0f : _rain * (HYDRO_ITER==0 ? 1.0f : 0.01f);
					
					if (rain < 0) {
						printf("rain = %f\n", rain);
					}
					*/
					
					int p = grid.At(_x, _y);
					
					grid.water[p] += 0.0001f * (HYDRO_ITER < 100 ? 500 : 1);//rain;
					
					SUM += grid.ground[p].layers[0] + grid.sediment[p];
				}
			}
			SUM_MATERIAL = SUM;
		}
		
		if (true) {
			grid.water[grid.At(1, 1)] += 0.1;
			grid.water[grid.At(1, 1)] += 0.1 * pow(sin(HYDRO_ITER/15.0f) + 1.0f, 2);
			for (int _y = 13; _y < 18; ++_y) {
				for (int _x = 498; _x < 503; ++_x) {
					grid.water[grid.At(15, 500)] += 0.1 * pow(sin(HYDRO_ITER/80.0f), 4);
				}
			}
		}
	}
	
	grid.FullCycle();

	if (updateHeights) {
		const auto now = std::chrono::steady_clock::now();
		auto dur = (now - lastUpdateHeightsTexture);
		if (dur > std::chrono::milliseconds(240)) {
			lastUpdateHeightsTexture = now;
			updateHeights = false;
			for (int _y = 0; _y < height; ++_y) {
				for (int _x = 0; _x < width; ++_x) {
// 						const int i = _x + _y * width;
					const int i = _x * height + _y;
					const int t = grid.At(_x, _y);
					float h = grid.ground[t].layers[0];// + grid.water[t] + grid.sediment[t];
// 					float h = t->sediment;
					h /= HYDRO_EROSION_Y_SCALE;
					
					if (h > -10000 && h < 50000) {
					} else {
						h = 0;
					}
					
					vertHeights[i].h = h;
					vertHeights[i].w = grid.water[t] / HYDRO_EROSION_Y_SCALE;
				}
			}
		}
	}

	if (updateWaterHeights) {
		const auto now = std::chrono::steady_clock::now();
		auto dur = (now - lastUpdateHeightsTexture);
		if (dur > std::chrono::milliseconds(240)) {
			lastUpdateHeightsTexture = now;
// 				updateWaterHeights = false;
			for (int _y = 0; _y < height; ++_y) {
				for (int _x = 0; _x < width; ++_x) {
// 						const int i = _x + _y * width;
					const int i = _x * height + _y;
					const int t = grid.At(_x, _y);
					float h = grid.ground[t].layers[0] + grid.ground[t].layers[1];//grid.water[t];// + grid.ground[t].Total() + grid.sediment[t];
// 						float h = grid.flux[t].fluxArray[0];// + grid.ground[t].Total() + grid.sediment[t];
					h /= HYDRO_EROSION_Y_SCALE;
					
					if (h > -10000 && h < 50000) {
					} else {
						h = 0;
					}
					
					vertHeights[i].h = h;
					vertHeights[i].w = grid.water[t] / HYDRO_EROSION_Y_SCALE;
				}
			}
		}
	}

	const auto b = std::chrono::steady_clock::now();
	hydroErosionDuration =
		std::chrono::nanoseconds(b - a).count() / 1'000'000.0f;
	if (HYDRO_ITER < 10) {
		averageHydroIterationDuration = hydroErosionDuration;
	} else if (HYDRO_ITER < 50) {
		averageHydroIterationDuration = std::min(hydroErosionDuration,
				averageHydroIterationDuration);
	} else {
		averageHydroIterationDuration =
			0.98f * averageHydroIterationDuration
			+
			0.02f * hydroErosionDuration;
	}
	
	HYDRO_ITER = HYDRO_ITER + 1;
}
