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

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz) {
	long long test_n = matrix_sz * matrix_sz;

	uint64_t *cr;
	cr = new uint64_t[test_n * branch_sz];
	for (int bid = 0; bid < branch_sz; bid++)
		for (int i = 0; i < test_n; i++) {
			cr[bid*test_n + i] = bid * matrix_sz;
		}

	auto start = clock_start();

	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);

	IntFp *mat_a = new IntFp[test_n * branch_sz];
	IntFp *mat_b = new IntFp[test_n * branch_sz];
	IntFp *mat_c = new IntFp[test_n * branch_sz];
	
	for (int bid = 0; bid < branch_sz; bid++)
		for(int i = 0; i < test_n; ++i) {
			mat_a[bid*test_n + i] = IntFp(1, ALICE);
			mat_b[bid*test_n + i] = IntFp(1, ALICE);
			mat_c[i] = IntFp((uint64_t)0, PUBLIC);
		}

	for (int bid = 0; bid < branch_sz; bid++)
		for(int i = 0; i < matrix_sz; ++i) {
			for(int j = 0; j < matrix_sz; ++j) {
				for(int k = 0; k < matrix_sz; ++k) {
					IntFp tmp = mat_a[bid*test_n+i*matrix_sz+j] * mat_b[bid*test_n+j*matrix_sz+k];
					mat_c[bid*test_n+i*matrix_sz+k] = mat_c[bid*test_n+i*matrix_sz+k] + tmp;
				}
			}
		}

	batch_reveal_check(mat_c, cr, test_n*branch_sz);
	auto timeuse = time_from(start);
	finalize_zk_arith<BoolIO<NetIO>>();
	cout << matrix_sz << "\t" << timeuse << " us\t" << party << " " << endl;
	std::cout << std::endl;

	delete[] cr;
	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;


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
	if(argc < 3) {
		std::cout << "usage: bin/arith/matrix_mul_arith PARTY PORT DIMENSION" << std::endl;
		return -1;
	} else if (argc == 3) {
		num = 10;
		branch = 10;
	} else {
		num = atoi(argv[3]);
		branch = atoi(argv[4]);
	}
	

	test_circuit_zk(ios, party, num, branch);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
