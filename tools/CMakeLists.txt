# Tools using the QuickEd Library

add_custom_target(tools)

# Generate Dataset Tool
add_executable(generate_dataset ${CMAKE_CURRENT_SOURCE_DIR}/generate_dataset/generate_dataset.c)
target_link_libraries(generate_dataset m) # Link math library
target_sources(generate_dataset PRIVATE ${CMAKE_SOURCE_DIR}/quicked_utils/src/commons.c)
add_dependencies(tools generate_dataset)

# Align Benchmark Tool
add_executable(align_benchmark)
target_include_directories(align_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark)
target_link_libraries(align_benchmark PRIVATE m) # Link math library
target_link_libraries(align_benchmark PRIVATE quicked)

# Add OpenMP lib
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(align_benchmark PRIVATE OpenMP::OpenMP_CXX)
endif()
if(OpenMP_C_FOUND)
    target_link_libraries(align_benchmark PRIVATE OpenMP::OpenMP_C)
endif()

add_subdirectory(align_benchmark/external/edlib EXCLUDE_FROM_ALL)
target_link_libraries(align_benchmark PRIVATE edlib)

file(GLOB ALIGN_BENCH_SRCS 
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/edit/edit_bpm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/edit/edit_dp.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/benchmark/benchmark_edit.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/benchmark/benchmark_check.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/benchmark/benchmark_utils.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/align_benchmark_params.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/align_benchmark.c
    ${CMAKE_CURRENT_SOURCE_DIR}/align_benchmark/score_matrix.c)

file(GLOB UTILS_SRCS
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/commons.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/mm_allocator.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/vector.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/dna_text.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/cigar.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/profiler_counter.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/profiler_timer.c
    ${PROJECT_SOURCE_DIR}/quicked_utils/src/sequence_buffer.c
)
target_sources(align_benchmark PRIVATE ${ALIGN_BMRK_SRCS} ${UTILS_SRCS} ${ALIGN_BENCH_SRCS})
add_dependencies(tools align_benchmark) # TODO: Uncomment when ready