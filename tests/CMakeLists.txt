# Tests for the QuickEd Library
add_custom_target(tests)

file(GLOB FILE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/*.c") # Get all .c files in this directory
foreach (FILE ${FILE_LIST})
    cmake_path(GET FILE STEM TEST) # Remove .c extension

    add_executable(${TEST} ${FILE})
    target_link_libraries(${TEST} quicked)
    set_target_properties(${TEST} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/tests")

    # Add the test for each file
    add_test(NAME test_${TEST} COMMAND $<TARGET_FILE:${TEST}>)
    add_dependencies(tests ${TEST}) # 'tests' target depends from that executable
endforeach ()