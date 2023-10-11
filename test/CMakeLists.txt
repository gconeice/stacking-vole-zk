#Testing macro
macro (add_test_executable_with_lib folder _name)
	add_executable(test_${folder}_${_name} "${folder}/${_name}.cpp")
	target_link_libraries(test_${folder}_${_name} ${EMP-OT_LIBRARIES} emp-zk)
endmacro()

macro (add_test_case folder _name)
	add_test_executable_with_lib(${folder} ${_name})
  	add_test(NAME ${folder}_${_name} COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${folder}_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")
endmacro()

macro (add_test_case_with_run folder _name)
	add_test_executable_with_lib(${folder} ${_name})
	add_test(NAME ${folder}_${_name} COMMAND "./run" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${folder}_${_name}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/")

endmacro()

# Test Boolean Stack
# add_test_case_with_run(bool stack_mat_mul)
add_test_case_with_run(bool stack_mat_mul_RO)
# add_test_case_with_run(bool stack_sha256)
add_test_case_with_run(bool stack_sha256_RO)

# Test Single Stack
add_test_case_with_run(arith unstack_single_disj_matmul)
add_test_case_with_run(arith stack_single_disj_matmul)
# add_test_case_with_run(arith stack_single_disj_matmul_RO)

# Test Single Stack Online
# add_test_case_with_run(arith unstack_single_disj_matmul_online)
# add_test_case_with_run(arith stack_single_disj_matmul_online)
# add_test_case_with_run(arith stack_single_disj_matmul_online_RO)


# Test Batched Stack
add_test_case_with_run(arith unstack_batched_disj_matmul)
add_test_case_with_run(arith stack_multi_single_disj_matmul)
# add_test_case_with_run(arith stack_multi_single_disj_matmul_RO)
add_test_case_with_run(arith stack_batched_matmul_v1)
