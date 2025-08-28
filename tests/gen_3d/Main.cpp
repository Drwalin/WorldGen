#include <chrono>
#include <thread>
#include <atomic>

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

#include "../../include/worldgen/Noises.hpp"

int width = 500 * 4;
int height = 500 * 4;

float horizontalScale = 0.04f;
float verticalScale = 7;
float noiseHorizontalScale = 1.0f / 256.0f;

struct Vertex {
	float h;
	uint8_t rgba[4];
};
Vertex *verts;

void ThreadFunction();

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
	gl::VBO vbo(sizeof(float) + 4 * sizeof(uint8_t), gl::ARRAY_BUFFER,
				gl::DYNAMIC_DRAW);
	verts = new Vertex[width * height];
	vbo.Init(width * height);
	vbo.Generate(verts, width * height);

	std::thread(ThreadFunction).detach();

	// Initiate VAO with VBO attributes
	gl::VAO vao(gl::TRIANGLES);
	vao.Init();
	vao.SetAttribPointer(vbo, shader.GetAttributeLocation("height"), 1,
						 gl::FLOAT, false, 0);
	vao.SetAttribPointer(vbo, shader.GetAttributeLocation("color"), 4,
						 gl::UNSIGNED_BYTE, true, 4);
	vao.BindElementBuffer(ebo, gl::DataType::UNSIGNED_INT);

	// Get uniform locations
	int modelLoc = shader.GetUniformLocation("model");
	int viewLoc = shader.GetUniformLocation("view");
	int projLoc = shader.GetUniformLocation("projection");

	shader.Use();
	glUniform2iv(shader.GetUniformLocation("size"), 1,
				 std::vector<int32_t>{width, height}.data());
	GL_CHECK_PUSH_ERROR;
	// 	shader.SetVec2(shader.GetUniformLocation("size"),
	// (*(glm::vec2*)(std::vector<int32_t>{width, height}.data())));
	shader.SetVec3(shader.GetUniformLocation("scale"),
				   glm::vec3(horizontalScale, verticalScale, horizontalScale));

	// Init camera position
	camera.ProcessMouseMovement(175, 50);
	camera.position = glm::vec3(-25, 65, -25);

	auto beg = std::chrono::steady_clock::now();
	int frames = 0;
	float fps = 0;
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
		printf("\r    fps: %.2f         ", fps);
		fflush(stdout);

		DefaultIterationStart();

		{
			static int count = 0;
			++count;
			if (count % 2 == 0) {
				vbo.Update(verts, 0, sizeof(Vertex) * width * height);
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

		// Draw VAO
		vao.Draw();

		DefaultIterationEnd();
	}

	gl::openGL.Destroy();
	glfwTerminate();

	return 0;
}

glm::ivec3 ColorGradient(int x, int y)
{
	thread_local wg::SimplexNoise colorSimplex(12342312);

	float n = colorSimplex.Noise2(glm::vec3{x * 10, 00, y * 10});
	
	glm::vec3 c = {n,n,n};
	c = glm::vec3(n * 230.0f);
	c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
	return c;
}

void ThreadFunction()
{
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			int i = x + y * width;
			verts[i] = {0, {0, 0, 0, 1}};
		}
	}

	wg::SimplexNoise simplex(1234);

	int CHUNK_SIZE = 64;

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
					
					x += 100;
					y += 100;
					
					y += 500;
					
					x += (simplex.Fbm(glm::vec2(-x/53+100,  y/53-1000), 3, 0.5, 2.3, false, false, 1.0f)-0.5) * 10;
					y += (simplex.Fbm(glm::vec2(+x/53-1000, -y/53+100), 3, 0.5, 2.3, false, false, 1.0f)-0.5) * 10;

					v.x =
						simplex.Terrain(glm::vec2{x, y} * noiseHorizontalScale,
										horizontalScale) *
						verticalScale;
					
					maxH = std::max(maxH, v.x);

// 					v.x *= sqrt(v.x / verticalScale);
// 					v.x = sqrt(v.x * verticalScale);
					float h = v.x;
					glm::ivec3 c = ColorGradient(
						_x, _y);

					verts[i] = {h,
								{(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 1}};
				}
			}
		}
		printf("Max height: %.3f\n", maxH);
		threadsDone++;
	};

	int threadsCount = std::thread::hardware_concurrency() - 1;

	for (int i = 1; i < threadsCount; ++i) {
		threads.push_back(std::thread(threadCalcFunc));
		threads.back().detach();
	}

	threadCalcFunc();

	while (threadsDone.load() < threads.size() + 1) {
		std::this_thread::sleep_for(std::chrono::milliseconds(60));
	}
	printf("\r Done!                         \n");
}
