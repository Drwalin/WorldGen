#pragma once

#include <cinttypes>

#include <glm/common.hpp>
#include <glm/fwd.hpp>

namespace wg
{
namespace hx
{
float Hash(glm::ivec3 p, int seed = 0);

float Hash(glm::ivec2 p, int seed = 0);
glm::vec2 Hash2(glm::ivec2 p, int seed = 0);

uint64_t mxrmx(uint64_t x);
}
}
