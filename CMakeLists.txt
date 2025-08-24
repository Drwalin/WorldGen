cmake_minimum_required(VERSION 3.8)
project(WorldGen)
cmake_policy(SET CMP0115 OLD)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WORLDGEN_3D_RENDERER "Build OpenGLWrapper examples" OFF)

set(
	WORLDGEN_OPENGLWRAPPER_DIR
    "../OpenGLWrapper"
    CACHE
    STRING
    "OpenGLWrapper directory"
)

include_directories(./include/)
aux_source_directory(./include/worldgen/ header_files)
aux_source_directory(./src/ source_files)
aux_source_directory(./thirdparty/lodepng/ lodepng)

add_library(WorldGen STATIC
	${header_files}
	${source_files}
	${lodepng}
)

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

