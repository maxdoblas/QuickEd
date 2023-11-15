cmake_minimum_required(VERSION 3.20)

project(quicked C)
set(CMAKE_C_STANDARD 11)

enable_testing()

# Set build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Set compiler flags if Debug build
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Build type: ${CMAKE_BUILD_TYPE}")
    add_compile_options(-Wall -Wextra -Wpedantic)

    # # Extra-hardcore debug mode 🔥
    # add_compile_options(-Wno-unused-parameter -Wfloat-equal -Wshadow -Wpointer-arith
    #                     -Wcast-align -Wstrict-prototypes -Wwrite-strings -Waggregate-return
    #                     -Wcast-qual -Wswitch-default -Wswitch-enum -Wconversion)

    # Enable address sanitizer
    add_compile_options(-fsanitize=address -fsanitize=undefined -fsanitize=shift -fsanitize=alignment
                        -fsanitize=signed-integer-overflow -fsanitize=bool -fsanitize=enum
                        -fsanitize=pointer-compare -fsanitize=pointer-overflow -fsanitize=builtin
                        -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address -fsanitize=undefined -fsanitize=shift -fsanitize=alignment
                        -fsanitize=signed-integer-overflow -fsanitize=bool -fsanitize=enum
                        -fsanitize=pointer-compare -fsanitize=pointer-overflow -fsanitize=builtin)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -flto)
endif()

# Set architecture to local machine architecture
add_compile_options(-march=native)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(quicked)
add_subdirectory(examples)
add_subdirectory(tools)