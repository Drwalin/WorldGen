#include <glm/geometric.hpp>
#include <thread>

#include <glm/fwd.hpp>

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

int width = 512;
int height = 512;

float horizontalScale = 0.1f;
float verticalScale = 15;
float noiseHorizontalScale = 0.01f * 0.5;

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
	shader.Load("../tests/gen_3d/vertex.glsl", "",
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

	while (!glfwWindowShouldClose(gl::openGL.window)) {
		DefaultIterationStart();

		vbo.Update(verts, 0, sizeof(Vertex) * width * height);

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

glm::ivec3 Color(int x, int y, float h)
{
	glm::vec3 c = {
	(wg::Noise::FractalBrownianMotion({x*horizontalScale,00,y*horizontalScale}, 3)+2.0f).x * 63.0f,
	(wg::Noise::FractalBrownianMotion({x*horizontalScale,10,y*horizontalScale}, 3)+2.0f).x * 63.0f,
	(wg::Noise::FractalBrownianMotion({x*horizontalScale,20,y*horizontalScale}, 3)+2.0f).x * 63.0f};
	
	
	if (h < 0.3) {
		c = {c.x,0,0};
	} else if (h < 0.7) {
		c = {0,c.y,0};
	} else {
		c = {0,0,c.z};
	}
	
	return c;
}

glm::ivec3 ColorGradient(int x, int y, glm::vec3 g)
{
	glm::vec3 ga = glm::abs(g);
	ga.x = 0;
	float sl = glm::dot(ga, ga);
	glm::vec3 c;
	
	float n = ((wg::Noise::FractalBrownianMotion({x*horizontalScale*10,00,y*horizontalScale*10}, 3)+1.0f).x * 0.5f);
	
	if (sl > 0.3) {
		c = glm::vec3(((n/2.0f) + 0.5f) * 230.0f);
		c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
	} else {
		c = glm::vec3(((n/2.0f) + 0.5f) * 230.0f);
		c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
		c.x /= 9.0f;
		c.z /= 9.0f;
	}
	
	glm::vec2 d{g.y, g.z};
	glm::vec2 s{0.5, 0.2};
	s = glm::normalize(s);
	if (glm::dot(s, d) > 0.001) {
		c *= 0.7;
	}
	
// 	c = glm::vec3(n * 230.0f);
// 	c = glm::clamp(c, glm::vec3(0), glm::vec3(255));
	
	return c;
}

void ThreadFunction()
{
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			int i = x + y * width;
			verts[i] = {0, {0,0,0, 1}};
		}
	}
	
	for (int _x = 0; _x < width; ++_x) {
		printf("\r p: %5i           ", _x);
		fflush(stdout);
		for (int _y = 0; _y < height; ++_y) {
			int i = _x + _y * width;
			glm::vec3 v;
// 			v.x = wg::Noise::fbm(glm::vec2{x,y}*noiseHorizontalScale, 1, 5) * verticalScale;
// 			float v1 = wg::Noise::fbm(glm::vec2{(x+0.1),y}*noiseHorizontalScale, 1, 5) * verticalScale;
// 			float v2 = wg::Noise::fbm(glm::vec2{x,(y+0.1)}*noiseHorizontalScale, 1, 5) * verticalScale;
			
// 			v.x = wg::Noise::NoiseV(glm::vec2{x,y}*noiseHorizontalScale) * verticalScale;
// 			float v1 = wg::Noise::NoiseV(glm::vec2{(x+0.1),y}*noiseHorizontalScale) * verticalScale;
// 			float v2 = wg::Noise::NoiseV(glm::vec2{x,(y+0.1)}*noiseHorizontalScale) * verticalScale;
			
			float x = _x;
			float y = _y;
			
			x -= 200;
			y += 250;
			
// 			x += 115;
// 			y += 30;
// 			
// 			x += 127;
// 			y += 2;
			
			v.x = wg::Noise::Terrain(glm::vec2{x,y}*noiseHorizontalScale, horizontalScale) * verticalScale;
			float v1 = wg::Noise::Terrain(glm::vec2{(x+0.1),y}*noiseHorizontalScale, horizontalScale) * verticalScale;
			float v2 = wg::Noise::Terrain(glm::vec2{x,(y+0.1)}*noiseHorizontalScale, horizontalScale) * verticalScale;
			
			v.y = (v1 - v.x) / (0.1 * horizontalScale);
			v.z = (v2 - v.x) / (0.1 * horizontalScale);
// 			v *= verticalScale;
			float h = v.x;
// 			float h = wg::Noise::fbm({x*0.01f,0,y*0.01f}, 1, 5);
			glm::ivec3 c = ColorGradient(x, y, v);//glm::clamp(glm::vec3((h+1.0f)*63.0f*1.5f), glm::vec3(0), glm::vec3(255));
			
			verts[i] = {h, {(uint8_t)c.x, (uint8_t)c.y, (uint8_t)c.z, 1}};
		}
	}
	printf("\r Done!                         \n");
}
