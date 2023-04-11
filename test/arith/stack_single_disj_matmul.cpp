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
	long long mul_sz = matrix_sz * matrix_sz * matrix_sz;

	uint64_t *res;
	res = new uint64_t[test_n];
	for (int i = 0; i < test_n; i++) res[i] = 0;

	auto start = clock_start();

	setup_zk_arith<BoolIO<NetIO>>(ios, threads, party);


	/*
	IntFp x(12, ALICE);
	IntFp y(13, ALICE);
	IntFp z(156, ALICE);
	__uint128_t delta = ZKFpExec::zk_exec->get_delta();

	std::cout << "DELTA = " << LOW64(delta) << std::endl;

	if (party==BOB) {
		std::cout << add_mod(mult_mod(x.value, y.value), mult_mod(delta, z.value)) << std::endl;
	} else {
		std::cout << "X = " << LOW64(x.value) << ' ' << HIGH64(x.value) << std::endl;
		std::cout << add_mod(z.value, add_mod(PR-mult_mod(x.value, HIGH64(y.value)), PR-mult_mod(y.value, HIGH64(x.value)))) << std::endl;
		std::cout << mult_mod(x.value, y.value) << std::endl;
	}
	*/

	// Input
	IntFp *mat_a = new IntFp[test_n];
	IntFp *mat_b = new IntFp[test_n];
	for (int i = 0; i < test_n; i++) {
		mat_a[i] = IntFp(1, ALICE);
		mat_b[i] = IntFp(1, ALICE);
	}

	// Left, Right, Output of MUL gates
	IntFp *mul_le = new IntFp[mul_sz];
	IntFp *mul_ri = new IntFp[mul_sz];
	IntFp *mul_ou = new IntFp[mul_sz];
	for (int i = 0; i < mul_sz; i++) {
		mul_le[i] = IntFp(1, ALICE);
		mul_ri[i] = IntFp(1, ALICE);
		mul_ou[i] = mul_le[i] * mul_ri[i];
	}

	IntFp one = IntFp(1, PUBLIC);

	if (party == ALICE) {
		uint64_t xxx;
		ZKFpExec::zk_exec->recv_data(&xxx, sizeof(uint64_t));
		std::cout << "WOW " << xxx << std::endl;
	} else {
		uint64_t yyy = 12345;
		ZKFpExec::zk_exec->send_data(&yyy, sizeof(uint64_t));		
	}

	/*
	for (int bid = 0; bid < branch_sz; bid++) {
		IntFp *tmp_mat_c = new IntFp[test_n];
		for (int i = 0; i < test_n; i++) tmp_mat_c[i] = IntFp(add_mod(pr-(bid+1)*matrix_sz, 0), PUBLIC); // set the witness of branch i

		// Matrix Multipilication
		for(int i = 0; i < matrix_sz; ++i) {
			for(int j = 0; j < matrix_sz; ++j) {
				for(int k = 0; k < matrix_sz; ++k) {
					IntFp tmp = mat_a[i*matrix_sz+j] * mat_b[j*matrix_sz+k];
					tmp_mat_c[i*matrix_sz+k] = tmp_mat_c[i*matrix_sz+k] + tmp;
				}
			}
		}

		// Mux		
		for (int i = 0; i < matrix_sz; i++)
			for (int j = 0; j < matrix_sz; j++) {
				IntFp tmp = tmp_mat_c[i*matrix_sz+j] * select_vec[bid];
				mat_c[i*matrix_sz+j] = mat_c[i*matrix_sz+j] + tmp;
			}
	}
	*/

	//batch_reveal_check(mat_c, res, test_n);	
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
