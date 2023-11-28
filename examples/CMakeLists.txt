# Examples on how to use the QuickEd Library
add_custom_target(examples)

# Plain examples
file(GLOB FILE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/*.c") # Get all .c files in this directory
foreach (FILE ${FILE_LIST})
    cmake_path(GET FILE STEM EXAMPLE) # Remove .c extension

    add_executable(${EXAMPLE} ${FILE})
    target_link_libraries(${EXAMPLE} quicked)
    set_target_properties(${EXAMPLE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/examples")

    # Add a test for each example. Although this is not proper testing, at least we can check that the example runs
    add_test(NAME test_example_${EXAMPLE} COMMAND $<TARGET_FILE:${EXAMPLE}>)
    add_dependencies(examples ${EXAMPLE}) # 'examples' target depends from that executable
endforeach ()

# C++ Binding Examples
file(GLOB FILE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/bindings/*.cpp") # Get all *.cpp files in this directory
foreach (FILE ${FILE_LIST})
    cmake_path(GET FILE STEM EXAMPLE) # Remove .cpp extension

    add_executable(binding_${EXAMPLE}_cpp ${FILE})
    target_link_libraries(binding_${EXAMPLE}_cpp quickedcpp)
    set_target_properties(binding_${EXAMPLE}_cpp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/examples/bindings")
    target_include_directories(binding_${EXAMPLE}_cpp PRIVATE "${PROJECT_SOURCE_DIR}/bindings/cpp")

    # Add a test for each example. Although this is not proper testing, at least we can check that the example runs
    add_test(NAME test_binding_${EXAMPLE}_cpp COMMAND $<TARGET_FILE:binding_${EXAMPLE}_cpp>)
    add_dependencies(examples binding_${EXAMPLE}_cpp) # 'examples' target depends from that executable
endforeach ()