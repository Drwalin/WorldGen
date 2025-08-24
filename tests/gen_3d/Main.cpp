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

int width = 512;
int height = 512;

float horizontalScale = 1.0f;
float verticalScale = 1.0f;

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
	camera.position = glm::vec3(-75, 175, -75);

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

void ThreadFunction()
{
	for (int x = 0; x < width; ++x) {
		for (int y = 0; y < height; ++y) {
			int i = x + y * width;
			float h = (sin(x/10.0) + cos(y/10.0) + 2) * 0.25f;
			uint8_t c = (h+2)*63.0;
			verts[i] = {h, {c,c,c, 1}};
		}
	}
}
