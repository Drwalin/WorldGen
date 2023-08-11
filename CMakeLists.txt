cmake_minimum_required(VERSION 3.20)
project(world_gen VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(GLOB_RECURSE WORLD_GEN_SOURCES "src/*.cpp")
add_executable(world_gen ${WORLD_GEN_SOURCES})

add_custom_target(
	run ALL
	COMMAND world_gen
	DEPENDS world_gen
	COMMENT "Running world generation"
)
set_target_properties(run PROPERTIES EXCLUDE_FROM_ALL TRUE)

file(
	COPY "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json"  
	DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}"
)
