# QuickEd library
add_library(quicked STATIC)

# Include directories
target_include_directories(quicked PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include) # Library includes
target_include_directories(quicked PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # Library root
target_link_libraries(quicked PUBLIC m) # Link math library

# Sources
file(GLOB QUICKED_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
)
file(GLOB UTILS_SRCS
    "${PROJECT_SOURCE_DIR}/utils/src/commons.c"
    "${PROJECT_SOURCE_DIR}/utils/src/mm_allocator.c"
    "${PROJECT_SOURCE_DIR}/utils/src/profiler_timer.c"
    "${PROJECT_SOURCE_DIR}/utils/src/profiler_counter.c"
    "${PROJECT_SOURCE_DIR}/utils/src/vector.c"
    "${PROJECT_SOURCE_DIR}/utils/src/dna_text.c"
    "${PROJECT_SOURCE_DIR}/utils/src/cigar.c"
)
target_sources(quicked PRIVATE ${QUICKED_SRCS} ${UTILS_SRCS}) # Add sources to target