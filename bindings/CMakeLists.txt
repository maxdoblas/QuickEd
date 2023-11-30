# QuickEd library C++ binding
add_library(quickedcpp STATIC cpp/quicked.cpp)
target_include_directories(quickedcpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
target_link_libraries(quickedcpp PRIVATE quicked)

# QuickEd library Python binding
add_subdirectory(python/pybind11)
pybind11_add_module(pyquicked python/quicked.cpp)
target_include_directories(pyquicked PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp) # Python binding is based on C++ binding
target_link_libraries(pyquicked PRIVATE quickedcpp)