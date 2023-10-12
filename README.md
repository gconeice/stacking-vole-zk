# Disjunction and Batched Disjunction for VOLE-based Zero Knowledge

This is the artifact for the paper: **Batchman and Robin:
Batched and Non-batched Branching for Interactive ZK** to be presented on ACM CCS 2023.

Eprint link: https://eprint.iacr.org/2023/1257

Base
=====
We acknowledge that our protocols based on QuickSilver repo available at: https://github.com/emp-toolkit/emp-zk.
In particular, we fork the repo and develop based on it.
We also tweak some emp libraries.

The file `sha256.txt` is obtained from https://homes.esat.kuleuven.be/~nsmart/MPC.

Set up Environments Including EMP Libraries
=====

You can simply use `sudo bash setup.sh`. Or (step-by-step),

0. `mkdir setup && cd setup`
1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python[3] install.py --deps --tool --ot --zk`
    1. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    2. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).
3. `cd ..`
4. `sudo apt install -y emacs iperf iftop clang`

Install and Build
=====

You can simply use `bash install.sh`. Or,

1. `mkdir build && cd build && CC=clang CXX=clang++ cmake ../ && make && cp ../sha256.txt ./`

We tested above methods already on a vanilla Ubuntu 22.04 machine.

Test
=====
We have the following tests:
1. Boolean single disjunction:
   1. bool_stack_mat_mul_RO: Boolean matrix multiplications, RO version.
   2. bool_stack_mat_mul: Boolean matrix multiplications, Lemma 5.4 version.
   3. bool_stack_sha256_RO: Repeating SHA2, RO version.
   4. bool_stack_sha256: Repeating SHA2, Lemma 5.4 version.

2. Arithmetic single disjunction:
   1. arith_stack_single_disj_matmul: Arithmetic matrix multiplications, Lemma 5.4 version.
   2. arith_stack_single_disj_matmul_RO: Arithmetic matrix multiplications, RO version.
   3. arith_stack_single_disj_matmul_online: Arithmetic matrix multiplications, Lemma 5.4 version, online cost only.
   4. arith_stack_single_disj_matmul_online_RO: Arithmetic matrix multiplications, RO version, online cost only.

3. Arithmetic batched disjunction:
   1. arith_stack_batched_matmul_v1: Arithmetic matrix multiplications, RO version.

4. Arithmetic baseline QuickSilver:
   1. arith_unstack_single_disj_matmul: Arithmetic matrix multiplications, single disjunction, QuickSilver.
   2. arith_unstack_single_disj_matmul_online: Arithmetic matrix multiplications, single disjunction, QuickSilver, online cost only.
   3. arith_unstack_batched_disj_matmul: Arithmetic matrix multiplications, batched disjunction, QuickSilver.

5. Repeating single disjunction:
   1. arith_stack_multi_single_disj_matmul: Arithmetic matrix multiplications, repeating single disjunction, Lemma 5.4 version.
   2. arith_stack_multi_single_disj_matmul_RO: Arithmetic matrix multiplications, repeating single disjunction, RO version.

**Note that** the auto script will help you compile executable files (under `build/bin/`):
   1. test_bool_stack_mat_mul_RO
   2. test_bool_stack_sha256_RO
   3. test_arith_stack_single_disj_matmul
   4. test_arith_unstack_single_disj_matmul
   5. test_arith_stack_batched_matmul_v1
   6. test_arith_unstack_batched_disj_matmul
   7. test_arith_stack_multi_single_disj_matmul


To obtain other executable files, please modify the file `test/CMakeLists.txt`.
       

# How to Simulate Network Setting

**We use `tc` command to simulate the network setting.**

     DEV=lo
     
     sudo tc qdisc del dev $DEV root
     
     sudo tc qdisc add dev $DEV root handle 1: tbf rate 1Gbit burst 100000 limit 10000
     
     sudo tc qdisc add dev $DEV parent 1:1 handle 10: netem delay 2msec

Change DEV to the network card you need (e.g., ens5).
Note that `sudo tc qdisc del dev $DEV root` needs to be executed before resetting the network.
It is used to clean the `tc` setting.

Both P and V need to restrict the network, note that for 30ms latency, both parties should be set to `delay 15msec`.

**You can use `iperf` to test the network throughput. Namely:**

P: `iperf -s`

V: `iperf -c [ip addr]`

**You can use `ping` to test the network latency. Namely:**

V: `ping [ip addr]`

**How to test the Comm. in the paper?**

We tested it using linux command `iftop` e.g., on the P's machine, execute:

   iftop -i ens5 -f 'port 12345'
