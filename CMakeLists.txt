cmake_minimum_required(VERSION 3.20)

project(quicked LANGUAGES C CXX)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

include(CTest)
enable_testing()

# Set build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# # Print message if compiler is Clang
# if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
#     message(WARNING "Compiler: Clang")
# endif()

# Set compiler flags if Debug build
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Build type: ${CMAKE_BUILD_TYPE}")
    add_compile_options(-O0 -g -Wall -Wextra -Werror)

    # Ignore known non-concerning warnings
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:-Wno-override-init>")

    # Enable code coverage
    add_compile_options(--coverage)
    add_link_options(--coverage)

    # Extra debug flags
    add_compile_options(-Wfloat-equal -Wformat=2 -Wshift-overflow -Wdouble-promotion -Wunused-macros
                        -Wwrite-strings -Wcast-qual -Wswitch-default -Wswitch-enum)

    add_compile_options("$<$<C_COMPILER_ID:GNU>:-Wduplicated-branches;-Wduplicated-cond;-Wlogical-op>")

    add_compile_options("$<$<COMPILE_LANGUAGE:C>:-Wstrict-prototypes>")
    #add_compile_options(-Wconversion) # TODO: We would like to enable this, but it's too noisy for now

    # Enable address sanitizer
    add_compile_options(-fsanitize=address -fsanitize=undefined -fno-sanitize-recover -fno-omit-frame-pointer)
    add_compile_options(-fno-sanitize=shift) # Ignore shift sanitizer (false positives, as we do some 64-bit shifts)
    add_link_options(-fsanitize=address -fsanitize=undefined)

    add_compile_options(-fstack-protector) # Enable stack protection

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -flto -fPIC)
endif()

# Set architecture to local machine architecture
add_compile_options(-march=native)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/utils/include) # To implcitly #include .h inside utils/src

add_subdirectory(quicked)
add_subdirectory(bindings)
add_subdirectory(examples)
add_subdirectory(tools)