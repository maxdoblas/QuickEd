# Tools using the QuickEd Library

add_custom_target(tools)

# Generate Dataset Tool
add_executable(generate_dataset ${CMAKE_CURRENT_SOURCE_DIR}/generate_dataset/generate_dataset.c)
target_link_libraries(generate_dataset m) # Link math library
target_sources(generate_dataset PRIVATE ${CMAKE_SOURCE_DIR}/utils/src/commons.c)
add_dependencies(tools generate_dataset)

# Align Benchmark Tool
add_executable(align_benchmark EXCLUDE_FROM_ALL) # TODO: Remove EXCLUDE_FROM_ALL when ready
file (GLOB_RECURSE ALIGN_BMRK_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/*.c)
target_include_directories(align_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark)
target_sources(align_benchmark PRIVATE ${ALIGN_BMRK_SRCS})
# add_dependencies(tools align_benchmark) # TODO: Uncomment when ready