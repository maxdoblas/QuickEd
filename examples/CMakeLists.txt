# Examples on how to use the QuickEd Library
add_custom_target(examples)

file(GLOB FILE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/*.c") # Get all .c files in this directory
foreach (FILE ${FILE_LIST})
    cmake_path(GET FILE STEM EXAMPLE) # Remove .c extension

    add_executable(${EXAMPLE} ${FILE})
    target_link_libraries(${EXAMPLE} quicked)
    set_target_properties(${EXAMPLE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/examples")

    # Add a test for each example. Although this is not proper testing, at least we can check that the example runs
    add_test(NAME test_example_${EXAMPLE} COMMAND $<TARGET_FILE:${EXAMPLE}>)
    add_dependencies(examples ${EXAMPLE}) # 'examples' target depends from that executable
endforeach ()

# C++ binding example + test
add_executable(example_cpp ${CMAKE_CURRENT_SOURCE_DIR}/bindings/quickedcpp.cpp)
target_link_libraries(example_cpp quickedcpp)
set_target_properties(example_cpp PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/examples/bindings")
target_include_directories(example_cpp PRIVATE ${CMAKE_SOURCE_DIR}/bindings/cpp)

add_test(NAME test_example_cpp COMMAND $<TARGET_FILE:example_cpp>)
add_dependencies(examples ${EXAMPLE})