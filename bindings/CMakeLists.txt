# QuickEd library C++ binding
add_library(quickedcpp STATIC cpp/quicked.cpp)
target_include_directories(quickedcpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
target_link_libraries(quickedcpp PRIVATE quicked)