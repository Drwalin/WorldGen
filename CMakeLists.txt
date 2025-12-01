cmake_minimum_required(VERSION 3.8)
project(WorldGen)
cmake_policy(SET CMP0115 OLD)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WORLDGEN_3D_RENDERER "Build OpenGLWrapper examples" OFF)
option(USE_OPENGLWRAPPER "Enable OpenGLWrapper" OFF)

if(WORLDGEN_3D_RENDERER)
	set(USE_OPENGLWRAPPER ON)
endif()

set(
	WORLDGEN_OPENGLWRAPPER_DIR
    "../OpenGLWrapper"
    CACHE
    STRING
    "OpenGLWrapper directory"
)

include_directories(./include/)
aux_source_directory(./include/worldgen/ wg_header_files)
aux_source_directory(./src/ wg_source_files)
aux_source_directory(./thirdparty/lodepng/ lodepng)

add_library(WorldGen STATIC
	${wg_header_files}
	${wg_source_files}
	${lodepng}
	thirdparty/OpenSimplex2/OpenSimplex2F/OpenSimplex2F.cpp
	thirdparty/OpenSimplex2/OpenSimplex2F/OpenSimplex2F.hpp
)
target_link_libraries(WorldGen OpenGLWrapper)
if(USE_OPENGLWRAPPER)
	target_compile_definitions(WorldGen PUBLIC -DWORLDGEN_ENABLE_OPENGLWRAPPER)
endif()

if(WORLDGEN_3D_RENDERER)
	add_subdirectory(${WORLDGEN_OPENGLWRAPPER_DIR} ./build/OpenGLWrapper)
	aux_source_directory(./tests/gen_3d/ gen_3d_source)
	add_executable(gen_3d
		${gen_3d_source}
		${WORLDGEN_OPENGLWRAPPER_DIR}/samples/Camera.hpp
		${WORLDGEN_OPENGLWRAPPER_DIR}/samples/DefaultCameraAndOtherConfig.hpp
	)
	target_link_libraries(gen_3d WorldGen OpenGLWrapper)
endif()

add_executable(gen_image_ydro_erosion tests/TestHydroErosion.cpp)
target_link_libraries(gen_image_ydro_erosion WorldGen)

