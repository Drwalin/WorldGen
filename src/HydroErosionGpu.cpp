#include <cassert>

#include "../include/worldgen/HydroErosion.hpp"
#include "HydroErosionPure.h"

Grid::GPUCompute::~GPUCompute()
{
	delete shaderCalcOutFlux;
	delete shaderUpdateWaterLevelAndVelocity;
	delete shaderErosionAndDepositionCalculation;
	delete shaderErosionAndDepositionUpdate;
	delete shaderSedimentTransportation;
	delete shaderSedimentTransportationUpdate;
	delete shaderThermalErosionCalculation;
	delete shaderThermalErosionUpdate;
	delete shaderEvaporation;
	delete shaderEvaporationUpdate;
	delete shaderSmooth;
	delete shaderSmoothUpdate;
	delete shaderUpdateRainAndRiver;
	delete shaderUpdateHeightTexture;

	delete vboGround;
	delete vboVelocity;
	delete vboFlux;
	delete vboWaterSedimentTemp;
}

static const char *prefix = R"(

layout (local_size_x = 16, local_size_y = 16, local_size_z = 1) in;
void main() {
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x >= width || y >= height) {
		return;
	}
	)";
static const char *postfix = R"((x, y);
}
)";

void Grid::GPUCompute::Init(int w, int h, Grid *grid)
{
	this->grid = grid;
	this->width = w;
	this->height = h;

	shaderCalcOutFlux = new gl::Shader();
	shaderUpdateWaterLevelAndVelocity = new gl::Shader();
	shaderErosionAndDepositionCalculation = new gl::Shader();
	shaderErosionAndDepositionUpdate = new gl::Shader();
	shaderSedimentTransportation = new gl::Shader();
	shaderSedimentTransportationUpdate = new gl::Shader();
	shaderThermalErosionCalculation = new gl::Shader();
	shaderThermalErosionUpdate = new gl::Shader();
	shaderEvaporation = new gl::Shader();
	shaderEvaporationUpdate = new gl::Shader();
	shaderSmooth = new gl::Shader();
	shaderSmoothUpdate = new gl::Shader();
	shaderUpdateRainAndRiver = new gl::Shader();
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
		std::string("#version 430 core\n\n") +
		globalMacrosVariableBySize +
		gl::Shader::LoadFileUseIncludes("../src/HydroErosionPure.h").c_str() +
		"\n";

	struct Pair {
		gl::Shader *shader;
		const char *name;
	} pairs[] = {
		{shaderUpdateRainAndRiver, "UpdateRainAndRiver"},
		{shaderCalcOutFlux, "CalcOutFlux"},
		{shaderUpdateWaterLevelAndVelocity, "UpdateWaterLevelAndVelocity"},
		{shaderErosionAndDepositionCalculation,
		 "ErosionAndDepositionCalculation"},
		{shaderErosionAndDepositionUpdate, "ErosionAndDepositionUpdate"},
		{shaderSedimentTransportation, "SedimentTransportation"},
		{shaderSedimentTransportationUpdate, "SedimentTransportationUpdate"},
		{shaderThermalErosionCalculation, "ThermalErosionCalculation"},
		{shaderThermalErosionUpdate, "ThermalErosionUpdate"},
		{shaderEvaporation, "Evaporation"},
		{shaderEvaporationUpdate, "EvaporationUpdate"},
		{shaderSmooth, "Smooth"},
		{shaderSmoothUpdate, "SmoothUpdate"}};
	for (auto &it : pairs) {
		it.shader->Unuse();
		it.shader->Compile(baseCode + prefix + it.name + postfix);
		it.shader->Use();
		BindBuffers();
		SetUniforms(it.shader);
		it.shader->Unuse();
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

void Grid::GPUCompute::CallCalcOutFlux() { CallShader(shaderCalcOutFlux); }
void Grid::GPUCompute::CallUpdateWaterLevelAndVelocity()
{
	CallShader(shaderUpdateWaterLevelAndVelocity);
}
void Grid::GPUCompute::CallErosionAndDepositionCalculation()
{
	CallShader(shaderErosionAndDepositionCalculation);
}
void Grid::GPUCompute::CallErosionAndDepositionUpdate()
{
	CallShader(shaderErosionAndDepositionUpdate);
}
void Grid::GPUCompute::CallSedimentTransportation()
{
	CallShader(shaderSedimentTransportation);
}
void Grid::GPUCompute::CallSedimentTransportationUpdate()
{
	CallShader(shaderSedimentTransportationUpdate);
}

void Grid::GPUCompute::CallThermalErosionCalculation()
{
	CallShader(shaderThermalErosionCalculation);
}
void Grid::GPUCompute::CallThermalErosionUpdate()
{
	CallShader(shaderThermalErosionUpdate);
}

void Grid::GPUCompute::CallEvaporation() { CallShader(shaderEvaporation); }
void Grid::GPUCompute::CallEvaporationUpdate()
{
	CallShader(shaderEvaporationUpdate);
}

void Grid::GPUCompute::CallSmooth() { CallShader(shaderSmooth); }
void Grid::GPUCompute::CallSmoothUpdate() { CallShader(shaderSmoothUpdate); }

void Grid::GPUCompute::CallShader(gl::Shader *shader)
{
	shader->Use();
	SetUniforms(shader);
	shader->DispatchRoundGroupNumbers(width, height, 1);
	shader->Unuse();
}

void Grid::GPUCompute::CallUpdateRainAndRiver()
{
	CallShader(shaderUpdateRainAndRiver);
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

void Grid::GPUCompute::UpdateGround(GroundLayers *data)
{
	gl::VBO *vbo = vboGround;
	vbo->Update(data - PADDING, 0, grid->elementsStorage * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateVelocity(Velocity *data)
{
	gl::VBO *vbo = vboVelocity;
	vbo->Update(data - PADDING, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateFlux(Flux *data)
{
	gl::VBO *vbo = vboFlux;
	vbo->Update(data - PADDING, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateWater(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, 0, grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::UpdateSediment(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
void Grid::GPUCompute::UpdateTemp1(float *data)
{
	gl::VBO *vbo = vboWaterSedimentTemp;
	vbo->Update(data - PADDING, 2ll * grid->elementsStorage * sizeof(float), grid->elementsStorage * sizeof(float));
}
