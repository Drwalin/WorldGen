#include <cassert>

#include "../include/worldgen/HydroErosion.hpp"
#include "HydroErosionPure.h"

Grid::GPUCompute::~GPUCompute()
{
	delete vboGround;
	delete vboWater;
	delete vboSediment;
	delete vboTemp1;
	delete vboVelocity;
	delete vboFlux;

	delete vboRiverSources;
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

	vboGround = new gl::VBO(sizeof(ground[0]), gl::SHADER_STORAGE_BUFFER,
							gl::DYNAMIC_DRAW);
	vboWater = new gl::VBO(sizeof(water[0]), gl::SHADER_STORAGE_BUFFER,
						   gl::DYNAMIC_DRAW);
	vboSediment = new gl::VBO(sizeof(sediment[0]), gl::SHADER_STORAGE_BUFFER,
							  gl::DYNAMIC_DRAW);
	vboTemp1 = new gl::VBO(sizeof(temp1[0]), gl::SHADER_STORAGE_BUFFER,
						   gl::DYNAMIC_DRAW);
	vboVelocity = new gl::VBO(sizeof(velocity[0]), gl::SHADER_STORAGE_BUFFER,
							  gl::DYNAMIC_DRAW);
	vboFlux = new gl::VBO(sizeof(flux[0]), gl::SHADER_STORAGE_BUFFER,
						  gl::DYNAMIC_DRAW);

	vboRiverSources =
		new gl::VBO(12, gl::SHADER_STORAGE_BUFFER, gl::DYNAMIC_DRAW);

	vboGround->Init(width * height + 1);
	vboWater->Init(width * height + 1);
	vboSediment->Init(width * height + 1);
	vboTemp1->Init(width * height + 1);
	vboVelocity->Init(width * height + 1);
	vboFlux->Init(width * height + 1);
	vboRiverSources->Init(width * height + 1);

	const std::string baseCode =
		std::string("#version 430 core\n\n") +
		gl::Shader::LoadFileUseIncludes("../src/HydroErosionPure.h").c_str() + "\n";

	struct Pair {
		gl::Shader *shader;
		const char *name;
	} pairs[] = {
		{&shaderCalcOutFlux, "CalcOutFlux"},
		{&shaderUpdateWaterLevelAndVelocity, "UpdateWaterLevelAndVelocity"},
		{&shaderErosionAndDepositionCalculation,
		 "ErosionAndDepositionCalculation"},
		{&shaderErosionAndDepositionUpdate, "ErosionAndDepositionUpdate"},
		{&shaderSedimentTransportation, "SedimentTransportation"},
		{&shaderSedimentTransportationUpdate, "SedimentTransportationUpdate"},
		{&shaderThermalErosionCalculation, "ThermalErosionCalculation"},
		{&shaderThermalErosionUpdate, "ThermalErosionUpdate"},
		{&shaderEvaporation, "Evaporation"},
		{&shaderSmooth, "Smooth"},
		{&shaderSmoothUpdate, "SmoothUpdate"}};
	for (auto &it : pairs) {
		shaderEvaporation.Unuse();
		it.shader->Compile(baseCode + prefix + it.name + postfix);
		it.shader->Use();
		BindBuffers();
		SetUniforms(it.shader);
		shaderEvaporation.Unuse();
	}
	shaderEvaporation.Unuse();

	shaderUpdateHeightTexture.Compile(baseCode + R"(

layout (binding = 0, rg32f) uniform image2D heightsTex;

layout(binding = 7) BUFFER(BufferBlock7, vec3, riverSources);

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
	vec2 h = vec2(g.layers[0] + g.layers[1], w);
	
	imageStore(heightsTex, ivec2(x, y), vec4(h, 0, 0));
}
)");
	shaderUpdateHeightTexture.Use();
	BindBuffers();
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, vboRiverSources->GetIdGL());
	SetUniforms(&shaderUpdateHeightTexture);
	shaderUpdateHeightTexture.Unuse();

	shaderUpdateRainAndRiver.Compile(baseCode + R"(

layout(binding = 7) BUFFER(BufferBlock7, vec3, riverSources);

layout (local_size_x = 16, local_size_y = 16, local_size_z = 1) in;
void main() {
	int x = int(gl_GlobalInvocationID.x);
	int y = int(gl_GlobalInvocationID.y);
	if (x >= width || y >= height) {
		return;
	}
	int id = At(x, y);
	float w = water[id];
	water[id] = w + 0.01 * dt;
}
)");
	shaderUpdateRainAndRiver.Use();
	BindBuffers();
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, vboRiverSources->GetIdGL());
	SetUniforms(&shaderUpdateRainAndRiver);
	shaderUpdateRainAndRiver.Unuse();
}

void Grid::GPUCompute::BindBuffers()
{
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, vboGround->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, vboWater->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, vboSediment->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, vboTemp1->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, vboVelocity->GetIdGL());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, vboFlux->GetIdGL());
}

void Grid::GPUCompute::SetUniforms(gl::Shader *shader)
{
	shader->Use();
	shader->SetFloat(shader->GetUniformLocation("hardness"),
					 {grid->hardness[0], grid->hardness[1]});
	shader->SetFloat(
		shader->GetUniformLocation("tangentOfAngleOfRecluse"),
		{grid->tangentOfAngleOfRecluse[0], grid->tangentOfAngleOfRecluse[1]});
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

void Grid::GPUCompute::CallCalcOutFlux() { CallShader(&shaderCalcOutFlux); }
void Grid::GPUCompute::CallUpdateWaterLevelAndVelocity()
{
	CallShader(&shaderUpdateWaterLevelAndVelocity);
}
void Grid::GPUCompute::CallErosionAndDepositionCalculation()
{
	CallShader(&shaderErosionAndDepositionCalculation);
}
void Grid::GPUCompute::CallErosionAndDepositionUpdate()
{
	CallShader(&shaderErosionAndDepositionUpdate);
}
void Grid::GPUCompute::CallSedimentTransportation()
{
	CallShader(&shaderSedimentTransportation);
}
void Grid::GPUCompute::CallSedimentTransportationUpdate()
{
	CallShader(&shaderSedimentTransportationUpdate);
}
void Grid::GPUCompute::CallThermalErosionCalculation()
{
	CallShader(&shaderThermalErosionCalculation);
}
void Grid::GPUCompute::CallThermalErosionUpdate()
{
	CallShader(&shaderThermalErosionUpdate);
}
void Grid::GPUCompute::CallEvaporation() { CallShader(&shaderEvaporation); }
void Grid::GPUCompute::CallSmooth() { CallShader(&shaderSmooth); }
void Grid::GPUCompute::CallSmoothUpdate() { CallShader(&shaderSmoothUpdate); }

void Grid::GPUCompute::CallShader(gl::Shader *shader)
{
	shader->Use();
	shader->DispatchRoundGroupNumbers(width, height, 1);
	shader->Unuse();
}

void Grid::GPUCompute::CallUpdateRainAndRiver()
{
	gl::Shader *shader = &shaderUpdateRainAndRiver;
	shader->Use();

	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, vboRiverSources->GetIdGL());

	shader->DispatchRoundGroupNumbers(width, height, 1);
	shader->Unuse();
}

void Grid::GPUCompute::UpdateHeightsTexture(gl::Texture *tex)
{
	gl::Shader *shader = &shaderUpdateHeightTexture;
	shader->Use();

	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, vboRiverSources->GetIdGL());

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
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateWater(float *data)
{
	gl::VBO *vbo = vboWater;
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateSediment(float *data)
{
	gl::VBO *vbo = vboSediment;
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateTemp1(float *data)
{
	gl::VBO *vbo = vboTemp1;
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateVelocity(Velocity *data)
{
	gl::VBO *vbo = vboVelocity;
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateFlux(Flux *data)
{
	gl::VBO *vbo = vboFlux;
	vbo->Update(data, 0, vbo->GetVertexCount() * vbo->VertexSize());
}
void Grid::GPUCompute::UpdateRiverSources(glm::vec3 *data, int amount)
{
	gl::VBO *vbo = vboRiverSources;
	vbo->Update(data, 0,
				std::min<int>(vbo->GetVertexCount(), amount) *
					vbo->VertexSize());
}
