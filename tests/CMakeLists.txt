# Tests for the QuickEd Library
add_custom_target(tests)

# Compile the test harness
add_executable(quicked_harness quicked_harness.c)
target_link_libraries(quicked_harness quicked)
add_dependencies(tests quicked_harness)

# Add the tests
add_test(NAME test_empty COMMAND $<TARGET_FILE:quicked_harness> "" "")
set_property(TEST test_empty PROPERTY PASS_REGULAR_EXPRESSION "ERROR: Tried to align an empty sequence")

add_test(NAME test_nonDNA COMMAND $<TARGET_FILE:quicked_harness> "GATC" "GATO" 1)

add_test(NAME test_l1000_n10000_e10 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/random_test.sh 1000 10000 10 quicked_harness)
set_property(TEST test_l1000_n10000_e10 PROPERTY
    ENVIRONMENT "BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_test(NAME test_l10000_n1000_e10 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/random_test.sh 10000 1000 10 quicked_harness)
set_property(TEST test_l10000_n1000_e10 PROPERTY
    ENVIRONMENT "BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_test(NAME test_l100000_n100_e10 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/random_test.sh 100000 100 10 align_benchmark)
set_property(TEST test_l100000_n100_e10 PROPERTY
    ENVIRONMENT "BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_test(NAME test_l1000000_n10_e10 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/random_test.sh 1000000 10 10 align_benchmark)
set_property(TEST test_l1000000_n10_e10 PROPERTY FAIL_REGULAR_EXPRESSION "INACCURATE SCORE")
set_property(TEST test_l1000000_n10_e10 PROPERTY
    ENVIRONMENT "BIN_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_test(NAME test_MiniION_align_benchmark COMMAND $<TARGET_FILE:align_benchmark> -i ${CMAKE_CURRENT_SOURCE_DIR}/test_data/ONT.MiniION.1.seq -c "score" -v)
set_property(TEST test_l1000000_n10_e10 PROPERTY FAIL_REGULAR_EXPRESSION "INACCURATE SCORE")