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
// #include "../../../OpenGLWrapper/include/openglwrapper/Sync.hpp"
// #include "../../../OpenGLWrapper/include/openglwrapper/Texture.hpp"
// #include "../../../OpenGLWrapper/include/openglwrapper/BufferAccessor.hpp"

#include "../../include/worldgen/HydroErosion.hpp"
#include "../../include/worldgen/Noises.hpp"

int width = 512 + 64;
int height = width;

float uniScale = 0.01;

float SCALE = 10;

float HYDRO_EROSION_Y_SCALE = 1.0f;

float horizontalScale = 1.0f * uniScale;
float verticalScale = uniScale;//1200.0f * uniScale / SCALE;
float noiseHorizontalScale = 1.0f / 4000.0f * SCALE;
volatile int HYDRO_ITER = 0;
volatile float averageHydroIterationDuration = 0.0f;
volatile bool disableSimulation = false;

volatile long double SUM_MATERIAL = 0;

std::mt19937_64 mt(std::random_device{}());

struct Vertex {
	float h;
	uint8_t rgba[4];
	float w;
};
Vertex *verts;

void ThreadFunction();

volatile float hydroErosionDuration = 0;
volatile bool updateHeights = true;
volatile bool updateWaterHeights = false;

int main(int argc, char **argv)
{
	DefaultsSetup();

	// Load shader
	gl::Shader shader;
	shader.Load("../tests/gen_3d/vertex.glsl", "../tests/gen_3d/geometry.glsl",
				"../tests/gen_3d/fragment.glsl");

	// TODO: load from argv

	// Generate index data
	gl::VBO ebo(sizeof(uint32_t), gl::ELEMENT_ARRAY_BUFFER, gl::DYNAMIC_DRAW);
	ebo.Init();

	{
		std::vector<uint32_t> Ebo;
		Ebo.reserve(3 * 2 * (width - 1) * (height - 1));
		for (int y = 0; y + 1 < height; ++y) {
			for (int x = 0; x + 1 < width; ++x) {
				int id0 = x + y * width;
				int id1 = id0 + 1;
				int id2 = id0 + width;
				int id3 = id1 + width;
				Ebo.push_back(id0);
				Ebo.push_back(id1);
				Ebo.push_back(id2);
				Ebo.push_back(id2);
				Ebo.push_back(id1);
				Ebo.push_back(id3);
			}
		}
		// Generate VBO from vertex data
		ebo.Generate(Ebo.data(), Ebo.size());
	}

	// Generate vertex data
	gl::VBO vbo(2 * sizeof(float) + 4 * sizeof(uint8_t), gl::ARRAY_BUFFER,
				gl::DYNAMIC_DRAW);
	verts = (Vertex*)vbo.InitMapPersistent(nullptr, width*height, 
// 					gl::DYNAMIC_STORAGE_BIT |
					gl::MAP_READ_BIT |
					gl::MAP_WRITE_BIT
// 					gl::MAP_PERSISTENT_BIT
// 					gl::MAP_COHERENT_BIT |
// 					gl::CLIENT_STORAGE_BIT
			);
// 	verts = new Vertex[width * height];
// 	vbo.Init(width * height);
// 	vbo.Generate(verts, width * height);

	std::thread(ThreadFunction).detach();

	// Initiate VAO with VBO attributes
	gl::VAO vao(gl::TRIANGLES);
	vao.Init();
	vao.SetAttribPointer(vbo, shader.GetAttributeLocation("height"), 1,
						 gl::FLOAT, false, 0);
	vao.SetAttribPointer(vbo, shader.GetAttributeLocation("color"), 4,
						 gl::UNSIGNED_BYTE, true, 4);
	vao.SetAttribPointer(vbo, shader.GetAttributeLocation("water"), 1,
						 gl::FLOAT, false, 8);
	vao.BindElementBuffer(ebo, gl::DataType::UNSIGNED_INT);

	// Get uniform locations
	int modelLoc = shader.GetUniformLocation("model");
	int viewLoc = shader.GetUniformLocation("view");
	int projLoc = shader.GetUniformLocation("projection");
	int useWaterLoc = shader.GetUniformLocation("useWater");
	int gridWidthLoc = shader.GetUniformLocation("gridWidth");

	shader.Use();
	glUniform2iv(shader.GetUniformLocation("size"), 1,
				 std::vector<int32_t>{width, height}.data());
	GL_CHECK_PUSH_ERROR;
	// 	shader.SetVec2(shader.GetUniformLocation("size"),
	// (*(glm::vec2*)(std::vector<int32_t>{width, height}.data())));
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
		printf("\r    fps: %.2f         (%i) hydro: %6.2f ms  (%6.2f ms)     mat = %Lf                            ", fps, HYDRO_ITER,
			   hydroErosionDuration, averageHydroIterationDuration, SUM_MATERIAL);
		fflush(stdout);

		DefaultIterationStart();
		
		if (gl::openGL.IsKeyDown('U')) {
			updateHeights = true;
		}
		if (gl::openGL.IsKeyDown('I')) {
			updateWaterHeights = true;
		}
		
		if (gl::openGL.WasKeyPressed('O')) {
			disableRender = !disableRender;
		}
		if (gl::openGL.WasKeyPressed('P')) {
			disableSimulation = !disableSimulation;
		}

		{
// 			constexpr int div = 32;
// 			static int count = 0;
// 			const int n = width * height;
// 			++count;
			
// 			if (count % 32 == 0)
// 			{
// 				vbo.Update(verts, 0, n * sizeof(Vertex));
				
// 				int mod = count % div;
// 				int of = (n * mod) / div;
// 				int ofn = ((n+1) * mod) / div;
// 				if (ofn >  n) {
// 					ofn = n;
// 				}
// 				int e = ofn - of;
// 				
// 				int bof = of * sizeof(Vertex);
				
// 				if (mod < 5) {
// 					vbo.Update(verts + of, bof, e * sizeof(Vertex));
// 				}
// 				gl::Finish();
// 			}
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

		// Draw VAO
		if (disableRender == false) {
			if (width > 1600) {
				glDepthMask(true);
				shader.SetInt(useWaterLoc, 1);
				vao.Draw();
			} else {
				shader.SetInt(useWaterLoc, 0);
				glDepthMask(true);
				vao.Draw();
				shader.SetInt(useWaterLoc, 1);
				glDepthMask(false);
				vao.Draw();
				glDepthMask(true);
			}
		}

		DefaultIterationEnd();
		
		gl::openGL.PrintErrors();
		
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}

	gl::openGL.Destroy();
	glfwTerminate();

	return 0;
}

glm::ivec3 ColorGradient(int x, int y)
{
	thread_local wg::SimplexNoise colorSimplex(12342312);

	float n = colorSimplex.Noise2(glm::vec3{x * 10, 00, y * 10});

	glm::vec3 c = {n, n, n};
	c = glm::vec3(n * 230.0f);
	c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
	return c;
}

void HydroErosionIteration();
Grid grid;

void ThreadFunction()
{
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			int i = x + y * width;
			verts[i] = {0, {0, 0, 0, 1}, 0};
		}
	}

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
					float x = _x;
					float y = _y;
					const int i = _x + _y * width;
					glm::vec3 v;

					// 					x += 100 + 500;
					// 					y += 100 + 500 + 300;
					//
					y += 750.f / SCALE;

					x += (simplex.Fbm(glm::vec2(-x / 53 + 100, y / 53 - 1000),
									  3, 0.5, 2.3, false, false, 1.0f) -
						  0.5) *
						 10;
					y += (simplex.Fbm(glm::vec2(+x / 53 - 1000, -y / 53 + 100),
									  3, 0.5, 2.3, false, false, 1.0f) -
						  0.5) *
						 10;

					v.x =
						simplex.Terrain(glm::vec2{x, y} * noiseHorizontalScale,
										horizontalScale);
					// 					v.x =
					// 						simplex.Noise(glm::vec2{x, y} *
					// noiseHorizontalScale) * 						verticalScale;

					maxH = std::max(maxH, v.x);

					// 					v.x *= sqrt(v.x / verticalScale);
					// 					v.x = sqrt(v.x * verticalScale);
					float h = v.x * 600.0f;
					glm::ivec3 c = ColorGradient(_x, _y);

					verts[i] = {h,
								{(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 1}, 0};
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

	HydroErosionIteration();
}

void HydroErosionIteration()
{
	grid.Init(width, height);
	wg::SimplexNoise simplex(432127);

	{
		wg::SimplexNoise simplex(13222);
		for (int _y = 0; _y < height; ++_y) {
			for (int _x = 0; _x < width; ++_x) {
				const int i = _x + _y * width;
				int t = grid.At<false>(_x, _y);
				grid.ground[t][0] = verts[i].h * HYDRO_EROSION_Y_SCALE - 2.0f;
				grid.ground[t][1] = 2.0f;
			}
		}
	}

	for (HYDRO_ITER=0;; HYDRO_ITER = HYDRO_ITER + 1) {
		if (disableSimulation) {
			std::this_thread::sleep_for(std::chrono::milliseconds(150));
			HYDRO_ITER = HYDRO_ITER - 1;;
			continue;
		}
		
		const auto a = std::chrono::steady_clock::now();
		
		for (int i=0; i<width * pow(width, 0.2); ++i) {
			int id = (mt()%(width*height)) + 1;
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
					
					int p = grid.At<false>(_x, _y);
					
					grid.water[p] += 0.0001f * (HYDRO_ITER < 100 ? 500 : 1);//rain;
					
					SUM += grid.ground[p][0] + grid.sediment[p];
				}
			}
			SUM_MATERIAL = SUM;
		}
		
		if (true) {
			grid.water[grid.At<false>(1, 1)] += 0.1;
			grid.water[grid.At<false>(1, 1)] += 0.1 * pow(sin(HYDRO_ITER/15.0f) + 1.0f, 2);
			for (int _y = 13; _y < 18; ++_y) {
				for (int _x = 498; _x < 503; ++_x) {
					grid.water[grid.At<false>(15, 500)] += 0.1 * pow(sin(HYDRO_ITER/80.0f), 4);
				}
			}
		}
// 		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		
		grid.FullCycle();
		
// 		std::this_thread::sleep_for(std::chrono::milliseconds(1500));

		if (updateHeights) {
			updateHeights = false;
			for (int _y = 0; _y < height; ++_y) {
				for (int _x = 0; _x < width; ++_x) {
					const int i = _x + _y * width;
					const int t = grid.At<false>(_x, _y);
					float h = grid.ground[t][0];// + grid.water[t] + grid.sediment[i];
// 					float h = t->sediment;
					h /= HYDRO_EROSION_Y_SCALE;
					
					if (h > -10000 && h < 50000) {
					} else {
						h = 0;
					}
					
					glm::ivec3 c = ColorGradient(_x, _y);
					verts[i] = {h, {(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 1}, grid.water[t] / HYDRO_EROSION_Y_SCALE};
				}
			}
		}

		if (updateWaterHeights) {
			updateWaterHeights = false;
			for (int _y = 0; _y < height; ++_y) {
				for (int _x = 0; _x < width; ++_x) {
					const int i = _x + _y * width;
					const int t = grid.At<false>(_x, _y);
					float h = grid.ground[t].Total();//grid.water[t];// + grid.ground[t].Total() + grid.sediment[i];
// 					float h = grid.flux[t].fluxArray[0];// + grid.ground[t].Total() + grid.sediment[i];
					h /= HYDRO_EROSION_Y_SCALE;
					
					if (h > -10000 && h < 50000) {
					} else {
						h = 0;
					}
					
					glm::ivec3 c = ColorGradient(_x, _y);
					verts[i] = {h, {(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 1}, grid.water[t] / HYDRO_EROSION_Y_SCALE};
				}
			}
		}
		
// 		std::this_thread::sleep_for(std::chrono::milliseconds(150));

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
	}
}
