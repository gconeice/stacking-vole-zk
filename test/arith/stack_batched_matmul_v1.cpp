#include "emp-zk/emp-zk.h"
#include <iostream>
#include "emp-tool/emp-tool.h"
#if defined(__linux__)
	#include <sys/time.h>
	#include <sys/resource.h>
#elif defined(__APPLE__)
	#include <unistd.h>
	#include <sys/resource.h>
	#include <mach/mach.h>
#endif

using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

inline uint64_t calculate_hash(PRP &prp, uint64_t x) {
	block bk = makeBlock(0, x);
	prp.permute_block(&bk, 1);
	return LOW64(bk) % PR;
}

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz, int batch_sz) {
	long long test_n = matrix_sz * matrix_sz;
	long long mul_sz = matrix_sz * matrix_sz * matrix_sz;
    long long test_n_batch_sz = test_n * batch_sz;
    long long mul_sz_batch_sz = mul_sz * batch_sz;
    long long w_length = test_n * 2 + mul_sz * 3 + 1;
    long long w_length_batch_sz = w_length * batch_sz;
    long long w_length_branch_sz = w_length * branch_sz;
    long long check_length = mul_sz * 2 + test_n;

	auto start = clock_start();

	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);

	// Input
	IntFp *mat_a = new IntFp[test_n_batch_sz];
	IntFp *mat_b = new IntFp[test_n_batch_sz];
	for (int i = 0; i < test_n_batch_sz; i++) {
		mat_a[i] = IntFp(1, ALICE);
		mat_b[i] = IntFp(1, ALICE);
	}

	// Left, Right, Output of MUL gates
	IntFp *mul_le = new IntFp[mul_sz_batch_sz];
	IntFp *mul_ri = new IntFp[mul_sz_batch_sz];
	IntFp *mul_ou = new IntFp[mul_sz_batch_sz];
	for (int i = 0; i < mul_sz_batch_sz; i++) {
		mul_le[i] = IntFp(1, ALICE);
		mul_ri[i] = IntFp(1, ALICE);
		mul_ou[i] = mul_le[i] * mul_ri[i];
	}
	ZKFpExec::zk_exec->flush_and_proofs();

	// Unit for constant offsets
	IntFp one = IntFp(1, PUBLIC);
	uint64_t delta = ZKFpExec::zk_exec->get_delta();

    // Bob generates linear challenge r on the left-hand-side by a seed
    block left_r_seed;
    if (party == ALICE) {
		ZKFpExec::zk_exec->recv_data(&left_r_seed, sizeof(block));
    } else {
        PRG().random_block(&left_r_seed, 1);
        ZKFpExec::zk_exec->send_data(&left_r_seed, sizeof(block));
    }
    PRG prg_left_r(&left_r_seed);
    uint64_t *left_r = new uint64_t[check_length];
    for (int i = 0; i < check_length; i++) {
        block tmp;
        prg_left_r.random_block(&tmp, 1);
        left_r[i] = LOW64(tmp) % PR;
    }

    // Alice and Bob calulate the left_vec_a
    uint64_t *left_vec_a = new uint64_t[w_length_branch_sz];
    // for each branch
    for (int bid = 0; bid < branch_sz; bid++) {
        uint64_t head = w_length * bid;
        // columns of inputs mat_a
        for (int i = 0; i < matrix_sz; i++)
            for (int j = 0; j < matrix_sz; j++) {
                left_vec_a[head + i*matrix_sz + j] = 0;
                for (int k = 0; k < matrix_sz; k++) left_vec_a[head + i*matrix_sz + j] = add_mod(left_vec_a[head + i*matrix_sz + j], left_r[i*test_n + j*matrix_sz + k]);
            }
        head = head + test_n;
        // columns of inputs mat_b
        for (int i = 0; i < matrix_sz; i++)
            for (int j = 0; j < matrix_sz; j++) {
                left_vec_a[head + i*matrix_sz + j] = 0;
                for (int k = 0; k < matrix_sz; k++) left_vec_a[head + i*matrix_sz + j] = add_mod(left_vec_a[head + i*matrix_sz + j], left_r[mul_sz + k*test_n + i*matrix_sz + j]);
            }
        head = head + test_n;
        // columns of mul_le and mul_ri
        // le
        for (int i = 0; i < mul_sz; i++) left_vec_a[head + i] = PR - left_r[i];
        head = head + mul_sz;
        // ri
        for (int i = 0; i < mul_sz; i++) left_vec_a[head + i] = PR - left_r[mul_sz + i];
        head = head + mul_sz;
        // columns of mul_ou
        for (int i = 0; i < matrix_sz; i++)
            for (int j = 0; j < matrix_sz; j++)
                for (int k = 0; k < matrix_sz; k++) 
                    left_vec_a[head + i*test_n + j*matrix_sz + k] = PR - left_r[mul_sz*2 + i*matrix_sz + k];
        head = head + mul_sz;
        // last columns of constant 1 offset
        left_vec_a[head] = 0;
        for (int i = 0; i < test_n; i++) left_vec_a[head] = add_mod(left_vec_a[head], mult_mod((bid+1)*matrix_sz, left_r[mul_sz*2 + i]));
    }

    

    // Alice chooses active left vector to prove inner_product
    IntFp *left_v = new IntFp[w_length_batch_sz];
    for (int i = 0; i < w_length_batch_sz; i++) left_v[i] = IntFp(1, ALICE);

	finalize_zk_arith<BoolIO<NetIO>>();
	auto timeuse = time_from(start);	
	cout << matrix_sz << "\t" << timeuse << " us\t" << party << " " << endl;
	std::cout << std::endl;


#if defined(__linux__)
	struct rusage rusage;
	if (!getrusage(RUSAGE_SELF, &rusage))
		std::cout << "[Linux]Peak resident set size: " << (size_t)rusage.ru_maxrss << std::endl;
	else std::cout << "[Linux]Query RSS failed" << std::endl;
#elif defined(__APPLE__)
	struct mach_task_basic_info info;
	mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
	if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS)
		std::cout << "[Mac]Peak resident set size: " << (size_t)info.resident_size_max << std::endl;
	else std::cout << "[Mac]Query RSS failed" << std::endl;
#endif
}

int main(int argc, char** argv) {
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:"127.0.0.1",port+i), party==ALICE);

	std::cout << std::endl << "------------ circuit zero-knowledge proof test ------------" << std::endl << std::endl;;

	int num = 0;
	int branch = 0;
    int batch = 0;
	if(argc < 3) {
		std::cout << "usage: bin/arith/matrix_mul_arith PARTY PORT DIMENSION" << std::endl;
		return -1;
	} else if (argc == 3) {
		num = 10;
		branch = 10;
        batch = 10;
	} else {
		num = atoi(argv[3]);
		branch = atoi(argv[4]);
        batch = atoi(argv[5]);
	}
	

	test_circuit_zk(ios, party, num, branch, batch);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
