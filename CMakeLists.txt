cmake_minimum_required(VERSION 3.20)

project(quicked C)
set(CMAKE_C_STANDARD 11)

include(CTest)
enable_testing()

# Set build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Set compiler flags if Debug build
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Build type: ${CMAKE_BUILD_TYPE}")
    add_compile_options(-g -Wall -Wextra)
    add_compile_options(-Wno-unknown-pragmas) # Ignore pragma region errors (for VSCode, code regions)
    add_compile_options(-Wno-override-init -Wno-maybe-uninitialized) # Ignore known non-concerning warnings

    # Enable code coverage (Doesn't work yet)
    # add_compile_options(-fprofile-arcs -ftest-coverage --coverage)
    # add_link_options(--coverage)
    # link_libraries(gcov)

    # Extra debug flags
    add_compile_options(-Wfloat-equal -Wshadow -Wformat=2 -Wlogical-op -Wshift-overflow=2
                        -Wcast-align -Wstrict-prototypes -Wwrite-strings #-Wconversion
                        -Wcast-qual -Wswitch-default -Wswitch-enum -Wduplicated-cond)

    # Enable address sanitizer
    add_compile_options(-fsanitize=address -fsanitize=undefined -fno-sanitize-recover -fno-omit-frame-pointer)
    add_compile_options(-fno-sanitize=shift) # Ignore shift sanitizer (false positives, as we do some 64-bit shifts)
    add_link_options(-fsanitize=address -fsanitize=undefined)

    add_compile_options(-fstack-protector) # Enable stack protection

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -flto)
endif()

# Set architecture to local machine architecture
add_compile_options(-march=native)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/utils/include) # To implcitly #include .h inside utils/src

add_subdirectory(quicked)
add_subdirectory(examples)
add_subdirectory(tools)