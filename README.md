# Disjunction and Batched Disjunction for VOLE-based Zero Knowledge

Base
=====
We acknowledge that our protocols based on QuickSilver repo available at: https://github.com/emp-toolkit/emp-zk.
In particular, we folk the repo and develop based on it.
We also tweat some emp libraries.
We will further clarify and obey the license in their repo if we open source.

Installation EMP libraries
=====
1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python[3] install.py --deps --tool --ot --zk`
    1. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    2. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).

Build
=====
1. `mkdir build && cd build && cmake ../ && make`

Test
=====
We have the following tests:
1. Boolean single disjunction:
   1. bool_stack_mat_mul_RO: Boolean matrix multiplications, RO version.
   2. bool_stack_mat_mul: Boolean matrix multiplications, Lemma 5.4 version.
   3. bool_stack_sha256_RO: Repeating SHA2, RO version.
   4. bool_stack_sah256: Repeating SHA2, Lemma 5.4 version.

2. Arithmetic single disjunction:
   1. arith_stack_single_disj_matmul: Arithmetic matrix multiplications, Lemma 5.4 version.
   2. arith_stack_single_disj_matmul_RO: Arithmetic matrix multiplications, RO version.
   3. arith_stack_single_disj_matmul_online: Arithmetic matrix multiplications, Lemma 5.4 version, online cost only.
   4. arith_stack_single_disj_matmul_online_RO: Arithmetic matrix multiplications, RO version, online cost only.

3. Arithmetic batched disjunction:
   1. arith_stack_batched_matmul_v1: Arithmetic matrix multiplications, RO version.

4. Arithmetic baseline QuickSilver:
   1. arith_unstack_single_disj_matmul: Arithmetic matrix multiplications, single disjunction, QuickSilver.
   2. arith_unstack_single_disj_matmul: Arithmetic matrix multiplications, single disjunction, QuickSilver, online cost only.
   3. arith_unstack_batched_disj_matmul: Arithmetic matrix multiplications, batched disjunction, QuickSilver.

5. Repeating single disjunction:
   1. arith_stack_multi_single_disj_matmul: Arithmetic matrix multiplications, repeating single disjunction, Lemma 5.4 version.
   2. arith_stack_multi_single_disj_matmul_RO: Arithmetic matrix multiplications, repeating single disjunction, RO version.

