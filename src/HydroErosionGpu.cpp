#include "../include/worldgen/HydroErosion.hpp"

#include "HydroErosionPure.h"
#include "../../OpenGLWrapper/include/openglwrapper/OpenGL.hpp"

Grid::GPUCompute::~GPUCompute()
{
	delete shaderUpdateHeightTexture;

	delete vboGround;
	delete vboVelocity;
	delete vboFlux;
	delete vboWaterSedimentTemp;
}

static const char *postfix = R"((x, y);
}
)";

void Grid::GPUCompute::Init(int w, int h, Grid *grid, int workGroupSizeDim)
{
	this->grid = grid;
	this->width = w;
	this->height = h;

	
	shaderUpdateHeightTexture = new gl::Shader();

	vboGround = new gl::VBO(sizeof(ground[0]), gl::SHADER_STORAGE_BUFFER,
							gl::DYNAMIC_DRAW);
	vboVelocity = new gl::VBO(sizeof(velocity[0]), gl::SHADER_STORAGE_BUFFER,
							  gl::DYNAMIC_DRAW);
	vboFlux = new gl::VBO(sizeof(flux[0]), gl::SHADER_STORAGE_BUFFER,
						  gl::DYNAMIC_DRAW);
	vboWaterSedimentTemp = new gl::VBO(sizeof(water_sediment_temp[0]), gl::SHADER_STORAGE_BUFFER,
						   gl::DYNAMIC_DRAW);

	vboGround->Init(grid->elementsStorage);
	vboVelocity->Init(grid->elementsStorage);
	vboFlux->Init(grid->elementsStorage);
	vboWaterSedimentTemp->Init(grid->elementsStorage);
	
	char globalMacrosVariableBySize[1024];
	snprintf(globalMacrosVariableBySize, sizeof(globalMacrosVariableBySize),
			"#define WIDTH %i\n"
			"#define HEIGHT %i\n"
			"#define ALIGNEMENT %i\n"
			"#define PADDING %i\n"
			"#define ELEMENTS %i\n"
			"#define ELEMENTS_STORAGE %i\n",
			width,
			height,
			ALIGNEMENT,
			PADDING,
			grid->elements,
			grid->elementsStorage);

	const std::string baseCode =
		std::string("#version 450 core\n\n") +
		globalMacrosVariableBySize +
		gl::Shader::LoadFileUseIncludes("../src/HydroErosionPure.h").c_str() +
		"\n";
	
	char prefix[4096];
	snprintf(prefix, sizeof(prefix)-1, R"(

layout (local_size_x = %i, local_size_y = %i, local_size_z = 1) in;
void main() {
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x >= width || y >= height) {
		return;
	}
	)", workGroupSizeDim, workGroupSizeDim);

	for (auto &st : grid->stages) {
		for (auto &it : st) {
			it.shader = new gl::Shader();
			it.shader->Unuse();
			it.shader->Compile(baseCode + prefix + it.functionName + postfix);
			it.shader->Use();
			BindBuffers();
			SetUniforms(it.shader);
			it.shader->Unuse();
			it.functionGpu = [this](gl::Shader *shader) {
				CallShader(shader);
			};
		}
	}

	shaderUpdateHeightTexture->Compile(baseCode + R"(

layout (binding = 0, rg32f) uniform image2D heightsTex;

layout (local_size_x = 16, local_size_y = 16, local_size_z = 1) in;
void main() {
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x >= width || y >= height) {
		return;
	}
	int id = At(x, y);
	GroundLayers g = ground[id];
	float w = water[id];
	float sed = sediment[id];
	vec2 h = vec2(g.layers[0] + g.layers[1] + sed, w);
// 	h.x = CalcSedimentCapacity(id, x, y) * 100.0;
// 	h.x = sed;
// 	if (x == 0) {
// 		h.x = y;
// 	} else {
// 		if (h.x < 0) {
// 			h.x *= 100.0;
// 		}
// 	}
	
	imageStore(heightsTex, ivec2(x, y), vec4(h, 0, 0));
}
)");
	shaderUpdateHeightTexture->Use();
	BindBuffers();
	SetUniforms(shaderUpdateHeightTexture);
	shaderUpdateHeightTexture->Unuse();
}

void Grid::GPUCompute::BindBuffers()
{
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, vboGround->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, vboVelocity->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, vboFlux->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, vboWaterSedimentTemp->GetIdGL());
}

void Grid::GPUCompute::SetUniforms(gl::Shader *shader)
{
	shader->Use();
	if (grid->iteration < 5) {
		shader->SetFloat(shader->GetUniformLocation("hardness"),
						 {grid->hardness[0], grid->hardness[1]});
		shader->SetFloat(shader->GetUniformLocation("tangentOfAngleOfRecluse"),
						 {grid->tangentOfAngleOfRecluse[0],
						  grid->tangentOfAngleOfRecluse[1]});
		shader->SetInt(shader->GetUniformLocation("width"), width);
		shader->SetInt(shader->GetUniformLocation("height"), height);
		shader->SetFloat(shader->GetUniformLocation("dt"), grid->dt);
		shader->SetFloat(shader->GetUniformLocation("A"), grid->A);
		shader->SetFloat(shader->GetUniformLocation("g"), grid->g);
		shader->SetFloat(shader->GetUniformLocation("l"), grid->l);
		shader->SetFloat(shader->GetUniformLocation("Kd"), grid->Kd);
		shader->SetFloat(shader->GetUniformLocation("Kc"), grid->Kc);
		shader->SetFloat(shader->GetUniformLocation("minimumSedimentCapacity"),
						 grid->minimumSedimentCapacity);
	}

	shader->SetInt(shader->GetUniformLocation("iteration"), grid->iteration);
	shader->Use();
	glUniform1i(shader->GetUniformLocation("iteration"), grid->iteration);
}

void Grid::GPUCompute::CallShader(gl::Shader *shader)
{
	gl::MemoryBarrier(gl::SHADER_STORAGE_BARRIER_BIT |
					  gl::ATOMIC_COUNTER_BARRIER_BIT);
	shader->Use();
	SetUniforms(shader);
	shader->DispatchRoundGroupNumbers(width, height, 1);
	shader->Unuse();
}

void Grid::GPUCompute::UpdateHeightsTexture(gl::Texture *tex)
{
	gl::Shader *shader = shaderUpdateHeightTexture;
	shader->Use();
	tex->Bind();
	GL_CHECK_PUSH_ERROR;
	glBindImageTexture(0, tex->GetIdGL(), 0, GL_FALSE, 0, GL_WRITE_ONLY,
					   gl::RG32F);
	GL_CHECK_PUSH_ERROR;

	shader->DispatchRoundGroupNumbers(width, height, 1);
	shader->Unuse();
}

void Grid::GPUCompute::UpdateAll()
{
	UpdateVelocity(grid->velocity);
	UpdateGround(grid->ground);
	UpdateFlux(grid->flux);
	
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(grid->water_sediment_temp, 0, vbo->GetBytes());
}
void Grid::GPUCompute::UpdateGround(const GroundLayers *data)
{
	gl::VBO *vbo = vboGround;
	vbo->Update(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::UpdateVelocity(const Velocity *data)
{
	gl::VBO *vbo = vboVelocity;
	vbo->Update(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::UpdateFlux(const Flux *data)
{
	gl::VBO *vbo = vboFlux;
	vbo->Update(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::UpdateWater(const float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, 0, grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::UpdateSediment(const float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::UpdateTemp1(const float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, 2ll * grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::UpdateTemp2(const float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, 3ll * grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}

void Grid::GPUCompute::FetchAll()
{
	FetchVelocity(grid->velocity);
	FetchGround(grid->ground);
	FetchFlux(grid->flux);
	
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Fetch(grid->water_sediment_temp, 0, vbo->GetBytes());
}
void Grid::GPUCompute::FetchGround(GroundLayers *data)
{
	gl::VBO *vbo = vboGround;
	vbo->Fetch(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::FetchVelocity(Velocity *data)
{
	gl::VBO *vbo = vboVelocity;
	vbo->Fetch(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::FetchFlux(Flux *data)
{
	gl::VBO *vbo = vboFlux;
	vbo->Fetch(data - PADDING, 0, vbo->GetBytes());
}
void Grid::GPUCompute::FetchWater(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Fetch(data - PADDING, 0, grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::FetchSediment(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Fetch(data - PADDING, grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::FetchTemp1(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Fetch(data - PADDING, 2ll * grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::FetchTemp2(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Fetch(data - PADDING, 3ll * grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
