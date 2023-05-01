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
// #define HVZK


inline uint64_t calculate_hash(PRP &prp, uint64_t x) {
	block bk = makeBlock(0, x);
	prp.permute_block(&bk, 1);
	return LOW64(bk) % PR;
}

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz, int batch_sz) {
	long long test_n = matrix_sz * matrix_sz;
	long long mul_sz = matrix_sz * matrix_sz * matrix_sz;

	/*
	uint64_t *res;
	res = new uint64_t[test_n];
	for (int i = 0; i < test_n; i++) res[i] = 0;
	// checking algo
	uint64_t *a;
	uint64_t *b;
	a = new uint64_t[test_n];
	b = new uint64_t[test_n];
	for (int i = 0; i < test_n; i++) a[i] = b[i] = i;
	for (int i = 0; i < matrix_sz; i++)
		for (int j = 0; j < matrix_sz; j++)
			for (int k = 0; k < matrix_sz; k++)
				res[i*matrix_sz+k] = add_mod(res[i*matrix_sz+k], mult_mod(a[i*matrix_sz+j], b[j*matrix_sz+k]));
	*/
	
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

	IntFp *mat_a = new IntFp[test_n];
	IntFp *mat_b = new IntFp[test_n];
	IntFp *mul_le = new IntFp[mul_sz];
	IntFp *mul_ri = new IntFp[mul_sz];
	IntFp *mul_ou = new IntFp[mul_sz];

    while (batch_sz--) {

        // Input
        for (int i = 0; i < test_n; i++) {
            mat_a[i] = IntFp(1, ALICE);
            mat_b[i] = IntFp(1, ALICE);
        }

        // Left, Right, Output of MUL gates
        for (int i = 0; i < mul_sz; i++) {
            mul_le[i] = IntFp(1, ALICE);
            mul_ri[i] = IntFp(1, ALICE);
            mul_ou[i] = mul_le[i] * mul_ri[i];
        }
        ZKFpExec::zk_exec->flush_and_proofs();

        // Unit for constant offsets
        IntFp one = IntFp(1, PUBLIC);
        //uint64_t chi; // random challenge
        //uint64_t delta = ZKFpExec::zk_exec->get_delta();
        //std::cout << "DELTA = " << delta << std::endl;

        IntFp com_cir[branch_sz];

        block s_seed; 
        if (party == ALICE) {
            ZKFpExec::zk_exec->recv_data(&s_seed, sizeof(block));
        } else {
            PRG().random_block(&s_seed, 1);
            ZKFpExec::zk_exec->send_data(&s_seed, sizeof(block));
        }
        PRG prg_s(&s_seed);           

        if (party == ALICE) {		

            uint64_t proofs[branch_sz];
            uint64_t values[branch_sz];
            // Bob needs to generate potential proofs for all branches
            for (int bid = 0; bid < branch_sz; bid++) {

                proofs[bid] = 0;
                values[bid] = 0;
                block tmptmp;
                prg_s.random_block(&tmptmp, 1);
    			uint64_t coeff = LOW64(tmptmp) % PR; 		

                // First 2*|C_MUL| rows to check that le/ri are formed correctly
                // Checking left
                for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
                    int i = mul_id / matrix_sz;
                    uint64_t tmp = add_mod(PR-LOW64(mul_le[mul_id].value), mat_a[i].value);
                    proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                    tmp = add_mod(PR-HIGH64(mul_le[mul_id].value), HIGH64(mat_a[i].value));
                    values[bid] = add_mod(values[bid], mult_mod(tmp, coeff));
                    prg_s.random_block(&tmptmp, 1);
    			    coeff = LOW64(tmptmp) % PR; // coeff as chi^k
                }
                // Checking right
                for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
                    int i = mul_id % test_n;
                    uint64_t tmp = add_mod(PR-LOW64(mul_ri[mul_id].value), mat_b[i].value);
                    proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                    tmp = add_mod(PR-HIGH64(mul_ri[mul_id].value), HIGH64(mat_b[i].value));
                    values[bid] = add_mod(values[bid], mult_mod(tmp, coeff));
                    prg_s.random_block(&tmptmp, 1);
    			    coeff = LOW64(tmptmp) % PR; // coeff as chi^k
                }	
                // Checking output: the active branch is the one that res[i][j] = matrix_sz
                for (int i = 0; i < matrix_sz; i++) 
                    for (int k = 0; k < matrix_sz; k++) {
                        uint64_t tmp = PR - LOW64(one.value);
                        uint64_t tmp2 = PR - HIGH64(one.value);
                        tmp = mult_mod(tmp, matrix_sz*(bid+1));
                        tmp2 = mult_mod(tmp2, matrix_sz*(bid+1));
                        for (int j = 0; j < matrix_sz; j++) {
                            tmp = add_mod(tmp, mul_ou[i*test_n + j*matrix_sz + k].value);
                            tmp2 = add_mod(tmp2, HIGH64(mul_ou[i*test_n + j*matrix_sz + k].value));
                        }
                        proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                        values[bid] = add_mod(values[bid], mult_mod(tmp2, coeff));
                        prg_s.random_block(&tmptmp, 1);
        			    coeff = LOW64(tmptmp) % PR; // coeff as chi^k			
                    }							
            }		

            for (int bid = 0; bid < branch_sz; bid++) {
                com_cir[bid].value = values[bid];
                com_cir[bid].value <<= 64;
                com_cir[bid].value ^= proofs[bid];
            }
        
        } else {
            uint64_t proofs[branch_sz];
            // Bob needs to generate potential proofs for all branches
            for (int bid = 0; bid < branch_sz; bid++) {

                proofs[bid] = 0;
                block tmptmp;
                prg_s.random_block(&tmptmp, 1);
    			uint64_t coeff = LOW64(tmptmp) % PR; 				

                // First 2*|C_MUL| rows to check that le/ri are formed correctly
                // Checking left
                for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
                    int i = mul_id / matrix_sz;
                    uint64_t tmp = add_mod(PR-LOW64(mul_le[mul_id].value), mat_a[i].value);
                    proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                    prg_s.random_block(&tmptmp, 1);
        			coeff = LOW64(tmptmp) % PR; // coeff as chi^k		
                }
                // Checking right
                for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
                    int i = mul_id % test_n;
                    uint64_t tmp = add_mod(PR-LOW64(mul_ri[mul_id].value), mat_b[i].value);
                    proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                    prg_s.random_block(&tmptmp, 1);
        			coeff = LOW64(tmptmp) % PR; // coeff as chi^k		
                }	
                // Checking output: the active branch is the one that res[i][j] = matrix_sz
                for (int i = 0; i < matrix_sz; i++) 
                    for (int k = 0; k < matrix_sz; k++) {
                        uint64_t tmp = PR - LOW64(one.value);
                        tmp = mult_mod(tmp, matrix_sz*(bid+1));
                        for (int j = 0; j < matrix_sz; j++) tmp = add_mod(tmp, mul_ou[i*test_n + j*matrix_sz + k].value);
                        proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
                        prg_s.random_block(&tmptmp, 1);
        			    coeff = LOW64(tmptmp) % PR; // coeff as chi^k						
                    }							
            }

            for (int bid = 0; bid < branch_sz; bid++) {
                com_cir[bid].value = proofs[bid];
            }

        }	

        IntFp final_res = IntFp(1, PUBLIC);
        for (int i = 0; i < branch_sz; i++) final_res = final_res * com_cir[i];
        batch_reveal_check_zero(&final_res, 1);
    }

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