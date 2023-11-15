# Examples on how to use the QuickEd Library
link_libraries(quicked)

add_custom_target(examples)

file(GLOB FILE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/*.c") # Get all .c files in this directory
foreach (FILE ${FILE_LIST})
    cmake_path(GET FILE STEM EXAMPLE) # Remove .c extension

    add_executable(${EXAMPLE} EXCLUDE_FROM_ALL ${FILE}) # We don't want to build it by default
    set_target_properties(${EXAMPLE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/examples")

    # Add a test for each example. Although this is not proper testing, at least we can check that the example runs
    add_test(NAME test_example_${EXAMPLE} COMMAND $<TARGET_FILE:${EXAMPLE}>)

    add_dependencies(examples ${EXAMPLE}) # 'examples' target depends from that executable
endforeach ()