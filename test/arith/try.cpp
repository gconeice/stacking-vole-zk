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

	IntFp *mat_a = new IntFp[test_n_batch_sz];
	IntFp *mat_b = new IntFp[test_n_batch_sz];
	IntFp *mul_le = new IntFp[mul_sz_batch_sz];
	IntFp *mul_ri = new IntFp[mul_sz_batch_sz];
	IntFp *mul_ou = new IntFp[mul_sz_batch_sz];

	auto start = clock_start();

	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);

	// Input
	for (int i = 0; i < test_n_batch_sz; i++) mat_a[i] = IntFp(1, ALICE);
	for (int i = 0; i < test_n_batch_sz; i++) mat_b[i] = IntFp(1, ALICE);

	// Left, Right, Output of MUL gates
	for (int i = 0; i < mul_sz_batch_sz; i++) mul_le[i] = IntFp(1, ALICE);
	for (int i = 0; i < mul_sz_batch_sz; i++) mul_ri[i] = IntFp(1, ALICE);
    for (int i = 0; i < mul_sz_batch_sz; i++)
		mul_ou[i] = mul_le[i] * mul_ri[i];
	ZKFpExec::zk_exec->flush_and_proofs();

	// Unit for constant offsets
	IntFp one = IntFp(1, PUBLIC);
	uint64_t delta = ZKFpExec::zk_exec->get_delta();

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
