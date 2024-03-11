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
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/commons.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/mm_allocator.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/profiler_timer.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/profiler_counter.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/vector.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/dna_text.c"
    "${PROJECT_SOURCE_DIR}/quicked_utils/src/cigar.c"
)
target_sources(quicked PRIVATE ${QUICKED_SRCS} ${UTILS_SRCS}) # Add sources to target

# Install
install(TARGETS quicked DESTINATION lib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/quicked.h DESTINATION include)
install(FILES
    ${PROJECT_SOURCE_DIR}/quicked_utils/include/commons.h
    ${PROJECT_SOURCE_DIR}/quicked_utils/include/mm_allocator.h
    ${PROJECT_SOURCE_DIR}/quicked_utils/include/profiler_timer.h
    ${PROJECT_SOURCE_DIR}/quicked_utils/include/profiler_counter.h
    ${PROJECT_SOURCE_DIR}/quicked_utils/include/vector.h
    DESTINATION include/quicked_utils/include
)