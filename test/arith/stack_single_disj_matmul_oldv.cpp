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

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz) {
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
	ZKFpExec::zk_exec->flush_and_proofs();

	// Unit for constant offsets
	IntFp one = IntFp(1, PUBLIC);
	uint64_t delta = ZKFpExec::zk_exec->get_delta();

	//std::cout << HIGH64(one.value) << ' ' << LOW64(one.value) << std::endl;

	if (party == ALICE) {
		uint64_t chi; // random challenge
		ZKFpExec::zk_exec->recv_data(&chi, sizeof(uint64_t));
		
		uint64_t active_proof = 0;
		uint64_t coeff = 1; // coeff as chi^k

		// Alice generates proof for the active branch
		// First 2*|C_MUL| rows to check that le/ri are formed correctly
		// Checking left
		for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
			int i = mul_id / matrix_sz;
			uint64_t tmp = add_mod(PR-LOW64(mul_le[mul_id].value), mat_a[i].value);
			active_proof = add_mod(active_proof, mult_mod(tmp, coeff));
			coeff = mult_mod(coeff, chi);
		}
		// Checking right
		for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
			int i = mul_id % test_n;
			uint64_t tmp = add_mod(PR-LOW64(mul_ri[mul_id].value), mat_b[i].value);
			active_proof = add_mod(active_proof, mult_mod(tmp, coeff));
			coeff = mult_mod(coeff, chi);
		}	
		// Checking output: the active branch is the one that res[i][j] = matrix_sz
		for (int i = 0; i < matrix_sz; i++) 
			for (int k = 0; k < matrix_sz; k++) {
				uint64_t tmp = PR - LOW64(one.value);
				tmp = mult_mod(tmp, matrix_sz);
				for (int j = 0; j < matrix_sz; j++) tmp = add_mod(tmp, mul_ou[i*test_n + j*matrix_sz + k].value);
				active_proof = add_mod(active_proof, mult_mod(tmp, coeff));
				coeff = mult_mod(coeff, chi);				
			}

		// setting fixed key aes
		block key;
		ZKFpExec::zk_exec->recv_data(&key, sizeof(block));
		PRP prp(key); //if no parameter is provided, a public, fixed key will be used
		
		// starting the proof
		// commit the proof of active branch
		active_proof = calculate_hash(prp, active_proof);
		IntFp com_active_proof(active_proof, ALICE);
		// get expected proofs of all branches from V
		uint64_t proofs[branch_sz];
		ZKFpExec::zk_exec->recv_data(proofs, branch_sz*sizeof(uint64_t));

		IntFp accumlate_term(1, PUBLIC);
		for (int bid = 0; bid < branch_sz; bid++) {
			IntFp term = com_active_proof + IntFp(proofs[bid], PUBLIC).negate();
			accumlate_term = accumlate_term * term;
		}

		// Alice commit her share of the final [0] 
		block com_final_key;
		PRG().random_block(&com_final_key,1);
		block com_final_pi = com_final_key ^ makeBlock(0, LOW64(accumlate_term.value));
		prp.permute_block(&com_final_pi, 1);
		ZKFpExec::zk_exec->send_data(&com_final_pi, sizeof(block));

#ifndef HVZK // only need for malicious BOB for selective failure attack
		ZKFpExec::zk_exec->recv_data(&delta, sizeof(uint64_t));
		// Alice recover Bob's keys
		// Input
		uint64_t *key_a = new uint64_t[test_n];
		uint64_t *key_b = new uint64_t[test_n];
		for (int i = 0; i < test_n; i++) {
			key_a[i] = add_mod(LOW64(mat_a[i].value), PR-delta);
			key_b[i] = add_mod(LOW64(mat_b[i].value), PR-delta);
		}

		// Left, Right, Output of MUL gates
		uint64_t *key_le = new uint64_t[mul_sz];
		uint64_t *key_ri = new uint64_t[mul_sz];
		uint64_t *key_ou = new uint64_t[mul_sz];
		for (int i = 0; i < mul_sz; i++) {
			key_le[i] = add_mod(LOW64(mul_le[i].value), PR-delta);
			key_ri[i] = add_mod(LOW64(mul_ri[i].value), PR-delta);
			key_ou[i] = add_mod(LOW64(mul_ou[i].value), PR-delta);
		}		

		// One for offsets
		uint64_t key_one = add_mod(LOW64(one.value), PR-delta);

		// Alice play Bob-in-her-head
		for (int bid = 0; bid < branch_sz; bid++) {

			uint64_t expected_proof = 0;
			coeff = 1; // coeff as chi^k			

			// First 2*|C_MUL| rows to check that le/ri are formed correctly
			// Checking left
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id / matrix_sz;
				uint64_t tmp = add_mod(PR-key_le[mul_id], key_a[i]);
				expected_proof = add_mod(expected_proof, mult_mod(tmp, coeff));
				coeff = mult_mod(coeff, chi);
			}
			// Checking right
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id % test_n;
				uint64_t tmp = add_mod(PR-key_ri[mul_id], key_b[i]);
				expected_proof = add_mod(expected_proof, mult_mod(tmp, coeff));
				coeff = mult_mod(coeff, chi);
			}	
			// Checking output: the active branch is the one that res[i][j] = matrix_sz
			for (int i = 0; i < matrix_sz; i++) 
				for (int k = 0; k < matrix_sz; k++) {
					uint64_t tmp = PR - key_one;
					tmp = mult_mod(tmp, matrix_sz*(bid+1));
					for (int j = 0; j < matrix_sz; j++) tmp = add_mod(tmp, key_ou[i*test_n + j*matrix_sz + k]);
					expected_proof = add_mod(expected_proof, mult_mod(tmp, coeff));
					coeff = mult_mod(coeff, chi);				
				}	

			// feed into fixed key aes to see whether it is the same as proofs[bid]						
			if (proofs[bid] != calculate_hash(prp, expected_proof)) error("Cheating Bob!");
		}
#endif

		// Alice opens her key
		ZKFpExec::zk_exec->send_data(&com_final_key, sizeof(block));

	} else {
		PRG tmpprg;
		uint64_t chi; // random challenge
		tmpprg.random_data(&chi, sizeof(uint64_t));
		chi = chi % PR;
		ZKFpExec::zk_exec->send_data(&chi, sizeof(uint64_t));	

		uint64_t proofs[branch_sz];

		// Bob needs to generate potential proofs for all branches
		for (int bid = 0; bid < branch_sz; bid++) {

			proofs[bid] = 0;
			uint64_t coeff = 1; // coeff as chi^k			

			// First 2*|C_MUL| rows to check that le/ri are formed correctly
			// Checking left
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id / matrix_sz;
				uint64_t tmp = add_mod(PR-LOW64(mul_le[mul_id].value), mat_a[i].value);
				proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
				coeff = mult_mod(coeff, chi);
			}
			// Checking right
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id % test_n;
				uint64_t tmp = add_mod(PR-LOW64(mul_ri[mul_id].value), mat_b[i].value);
				proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
				coeff = mult_mod(coeff, chi);
			}	
			// Checking output: the active branch is the one that res[i][j] = matrix_sz
			for (int i = 0; i < matrix_sz; i++) 
				for (int k = 0; k < matrix_sz; k++) {
					uint64_t tmp = PR - LOW64(one.value);
					tmp = mult_mod(tmp, matrix_sz*(bid+1));
					for (int j = 0; j < matrix_sz; j++) tmp = add_mod(tmp, mul_ou[i*test_n + j*matrix_sz + k].value);
					proofs[bid] = add_mod(proofs[bid], mult_mod(tmp, coeff));
					coeff = mult_mod(coeff, chi);				
				}							
		}

		// Set up fixed key aes hash function
		block key;
		PRG().random_block(&key,1);
		ZKFpExec::zk_exec->send_data(&key, sizeof(block));		
		PRP prp(key); //if no parameter is provided, a public, fixed key will be used

		// get Alice commitment on active branch
		IntFp com_active_proof(0, ALICE);
		// send Alice all expected proofs
		for (int i = 0; i < branch_sz; i++) proofs[i] = calculate_hash(prp, proofs[i]);
		ZKFpExec::zk_exec->send_data(proofs, branch_sz*sizeof(uint64_t));

		IntFp accumlate_term(1, PUBLIC);
		for (int bid = 0; bid < branch_sz; bid++) {
			IntFp term = com_active_proof + IntFp(proofs[bid], PUBLIC).negate();
			accumlate_term = accumlate_term * term;
		}

		// Get Alice commitment for the final proof
		block com_final_pi;
		ZKFpExec::zk_exec->recv_data(&com_final_pi, sizeof(block));

#ifndef HVZK // only need for malicious BOB for selective failure attack
		// Bob opens secret key \Delta
		ZKFpExec::zk_exec->send_data(&delta, sizeof(uint64_t));
#endif		

		// Alice open the proof
		block com_final_key;
		ZKFpExec::zk_exec->recv_data(&com_final_key, sizeof(block));
		com_final_key = com_final_key ^ makeBlock(0, LOW64(accumlate_term.value));
		prp.permute_block(&com_final_key, 1);

		if (HIGH64(com_final_pi)!=HIGH64(com_final_key)  || LOW64(com_final_pi)!=LOW64(com_final_key)) error("Invalid proof");

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
