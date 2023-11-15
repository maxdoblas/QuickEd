# QuickEd library
add_library(quicked STATIC)

# Include directories
target_include_directories(quicked PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include # Library includes
    ${PROJECT_SOURCE_DIR}/utils/include # To implcitly #include .h inside utils/src
)

target_include_directories(quicked PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # Library root

# Sources
file(GLOB QUICKED_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
)
file(GLOB UTILS_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/../utils/src/mm_allocator.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/../utils/src/vector.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/../utils/src/dna_text.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/../utils/src/cigar.c"
)
target_sources(quicked PRIVATE ${QUICKED_SRCS} ${UTILS_SRCS}) # Add sources to target