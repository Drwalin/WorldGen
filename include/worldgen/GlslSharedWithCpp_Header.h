#ifndef C_GLSL_SHARED_CODE_H
#define C_GLSL_SHARED_CODE_H

#if GL_core_profile
#define inline
#define IGNORE_UNUSED(VAR)
#define BUFFER(BLOCK, TYPE, VAR) buffer BLOCK { TYPE _padding##VAR[PADDING]; TYPE VAR[]; };
#define static
#else
%:include <glm/glm.hpp>
%:include <glm/vec2.hpp>
%:include <glm/vec3.hpp>
%:include <glm/vec4.hpp>
%:include <glm/mat2x2.hpp>
%:include <glm/mat3x3.hpp>
%:include <glm/mat4x4.hpp>
#define layout(...)
#define uniform
#define IGNORE_UNUSED(VAR) (void)VAR;
#define BUFFER(BLOCK, TYPE, VAR) static TYPE *VAR = nullptr;
#endif

#endif
