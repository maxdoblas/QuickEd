cmake_minimum_required(VERSION 3.20)

project(quicked LANGUAGES C CXX)
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

include(CTest)
enable_testing()

# Options
option(QUICKED_NONATIVE "Compile QuickEd for generic x86_64" OFF)
option(QUICKED_FORCESCALAR "Compile QuickEd without vector extensions (SSE,AVX)" OFF)

# Set build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

option(ASAN "Enable ASAN and UBSAN" OFF)

# Set compiler flags if Debug build
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Build type: ${CMAKE_BUILD_TYPE}")
    add_compile_options(-O0 -g -Wall -Wextra -Werror)

    # Ignore known non-concerning warnings
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:-Wno-override-init>")
    add_compile_options(-Wno-error=cpp) # #warning are not errors

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
    if (ASAN)
        add_compile_options(-fsanitize=address -fsanitize=undefined -fno-sanitize-recover -fno-omit-frame-pointer)
        add_compile_options(-fno-sanitize=shift) # Ignore shift sanitizer (false positives, as we do some 64-bit shifts)
        add_link_options(-fsanitize=address -fsanitize=undefined)
    endif()

    # Enable stack protection
    add_compile_options(-fstack-protector)

    # Add runtime checkers
    add_compile_definitions(_FORTIFY_SOURCE=1)

elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -flto=auto)
endif()

# Set global compiler flags
add_compile_options(-fPIC)
if (NOT QUICKED_NONATIVE)
    add_compile_options(-march=native)
else()
    add_compile_options(-march=x86-64)
    if (NOT QUICKED_FORCESCALAR)
        add_compile_options(-msse4.1 -mavx2)
    endif()
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/quicked_utils/include) # To implcitly #include .h inside quicked_utils/src

add_subdirectory(quicked)
add_subdirectory(bindings)
add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(tools)
