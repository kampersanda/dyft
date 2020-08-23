set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(sample_bin sample_bin.cpp)
target_link_libraries(sample_bin ${Boost_LIBRARIES})

add_executable(sample_int sample_int.cpp)
target_link_libraries(sample_int ${Boost_LIBRARIES})

add_executable(simhash simhash.cpp)
target_link_libraries(simhash ${Boost_LIBRARIES})

add_executable(precompute precompute.cpp)
target_link_libraries(precompute ${Boost_LIBRARIES})

add_executable(bvecs_to_bin bvecs_to_bin.cpp)
target_link_libraries(bvecs_to_bin ${Boost_LIBRARIES})

add_executable(gen_uniform gen_uniform.cpp)
target_link_libraries(gen_uniform ${Boost_LIBRARIES})


# # # # # # # # # # # # # # # # 
#   build_index_bin
# # # # # # # # # # # # # # # # 
add_executable(build_index_bin_dyft build_index_bin.cpp)
target_link_libraries(build_index_bin_dyft ${Boost_LIBRARIES})
set_target_properties(build_index_bin_dyft PROPERTIES COMPILE_DEFINITIONS "ALGO_DYFT")

add_executable(build_index_bin_hms1v build_index_bin.cpp)
target_link_libraries(build_index_bin_hms1v ${Boost_LIBRARIES})
set_target_properties(build_index_bin_hms1v PROPERTIES COMPILE_DEFINITIONS "ALGO_HMS1V")


# # # # # # # # # # # # # # # # 
#   build_index_int
# # # # # # # # # # # # # # # # 
add_executable(build_index_int_dyft build_index_int.cpp)
target_link_libraries(build_index_int_dyft ${Boost_LIBRARIES})
set_target_properties(build_index_int_dyft PROPERTIES COMPILE_DEFINITIONS "ALGO_DYFT")

add_executable(build_index_int_hms1dv build_index_int.cpp)
target_link_libraries(build_index_int_hms1dv ${Boost_LIBRARIES})
set_target_properties(build_index_int_hms1dv PROPERTIES COMPILE_DEFINITIONS "ALGO_HMS1DV")

add_executable(build_index_int_gv build_index_int.cpp)
target_link_libraries(build_index_int_gv ${Boost_LIBRARIES})
set_target_properties(build_index_int_gv PROPERTIES COMPILE_DEFINITIONS "ALGO_GV")


# # # # # # # # # # # # # # # # 
#   range_search_bin
# # # # # # # # # # # # # # # # 
add_executable(range_search_bin_ls range_search_bin.cpp)
target_link_libraries(range_search_bin_ls ${Boost_LIBRARIES})
set_target_properties(range_search_bin_ls PROPERTIES COMPILE_DEFINITIONS "ALGO_LS")

add_executable(range_search_bin_dyft range_search_bin.cpp)
target_link_libraries(range_search_bin_dyft ${Boost_LIBRARIES})
set_target_properties(range_search_bin_dyft PROPERTIES COMPILE_DEFINITIONS "ALGO_DYFT")

add_executable(range_search_bin_hms1v range_search_bin.cpp)
target_link_libraries(range_search_bin_hms1v ${Boost_LIBRARIES})
set_target_properties(range_search_bin_hms1v PROPERTIES COMPILE_DEFINITIONS "ALGO_HMS1V")


# # # # # # # # # # # # # # # # 
#   range_search_int
# # # # # # # # # # # # # # # # 
add_executable(range_search_int_ls range_search_int.cpp)
target_link_libraries(range_search_int_ls ${Boost_LIBRARIES})
set_target_properties(range_search_int_ls PROPERTIES COMPILE_DEFINITIONS "ALGO_LS")

add_executable(range_search_int_dyft range_search_int.cpp)
target_link_libraries(range_search_int_dyft ${Boost_LIBRARIES})
set_target_properties(range_search_int_dyft PROPERTIES COMPILE_DEFINITIONS "ALGO_DYFT")

add_executable(range_search_int_hms1dv range_search_int.cpp)
target_link_libraries(range_search_int_hms1dv ${Boost_LIBRARIES})
set_target_properties(range_search_int_hms1dv PROPERTIES COMPILE_DEFINITIONS "ALGO_HMS1DV")

add_executable(range_search_int_gv range_search_int.cpp)
target_link_libraries(range_search_int_gv ${Boost_LIBRARIES})
set_target_properties(range_search_int_gv PROPERTIES COMPILE_DEFINITIONS "ALGO_GV")
