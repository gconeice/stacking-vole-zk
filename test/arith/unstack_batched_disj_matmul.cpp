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

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz, int batch_sz) {
	long long test_n = matrix_sz * matrix_sz;

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

	IntFp *select_vec = new IntFp[branch_sz];
    IntFp *mat_a = new IntFp[test_n];
    IntFp *mat_b = new IntFp[test_n];
    IntFp *mat_c = new IntFp[test_n];

    for (int bc = 0; bc < batch_sz; bc++) {
        // Demux/mux vector
        select_vec[0] = IntFp(1, ALICE);
        for (int i = 1; i < branch_sz; i++) select_vec[i] = IntFp(0, ALICE);

        for (int i = 0; i < test_n; i++) {
            mat_a[i] = IntFp(1, ALICE);
            mat_b[i] = IntFp(1, ALICE);
            mat_c[i] = IntFp(0, PUBLIC);
        }

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
            delete[] tmp_mat_c;
        }

        batch_reveal_check(mat_c, res, test_n);	
    }
	finalize_zk_arith<BoolIO<NetIO>>();
	auto timeuse = time_from(start);	
	cout << matrix_sz << "\t" << timeuse << " us\t" << party << " " << endl;
	std::cout << std::endl;

	delete[] res;
	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;
	delete[] select_vec;

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

	if(argc < 7) {
		std::cout << "usage: a.out PARTY PORT IP DIMENSION #BRANCH #REPETITION" << std::endl;
		return -1;
	}
	
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:argv[3],port+i), party==ALICE);

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
        batch = 0;
	} else {
		num = atoi(argv[4]);
		branch = atoi(argv[5]);
        batch = atoi(argv[6]);
	}
	

	test_circuit_zk(ios, party, num, branch, batch);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
