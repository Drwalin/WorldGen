#if GL_core_profile
#define inline
#define IGNORE_UNUSED(VAR)
#define BUFFER(BLOCK, TYPE, VAR)                                               \
	buffer BLOCK { TYPE VAR[]; };
#define static
#else
#include <glm/glm.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#define layout(a)
#define uniform
#define IGNORE_UNUSED(VAR) (void)VAR;
#define BUFFER(BLOCK, TYPE, VAR) static TYPE *VAR = nullptr;
#endif
