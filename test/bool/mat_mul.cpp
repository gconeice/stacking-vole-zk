#include <emp-zk/emp-zk.h>
#include <iostream>
#include "emp-tool/emp-tool.h"
using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz) {

    uint64_t test_n = matrix_sz * matrix_sz;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
	auto start = clock_start();

	Bit *mat_a = new Bit[test_n];
	Bit *mat_b = new Bit[test_n];
	Bit *mat_c = new Bit[test_n];
	
	for(int i = 0; i < test_n; ++i) {
		mat_a[i] = Bit(true, ALICE);
		mat_b[i] = Bit(true, ALICE);
		mat_c[i] = Bit(false, PUBLIC);
	}

    std::cout << HIGH64(mat_c[0].bit) << ' ' << LOW64(mat_c[0].bit) << std::endl;

    block delta = get_bool_delta<BoolIO<NetIO>>(party);
    cout << getLSB(delta) << endl;
    std::cout << getLSB(mat_a[0].bit) << std::endl;

	for(int i = 0; i < matrix_sz; ++i) {
		for(int j = 0; j < matrix_sz; ++j) {
			for(int k = 0; k < matrix_sz; ++k) {
				Bit tmp = CircuitExecution::circ_exec->and_gate(mat_a[i*matrix_sz+j].bit, mat_b[j*matrix_sz+k].bit);                
                //mat_a[i*matrix_sz+j] & mat_b[j*matrix_sz+k];
				mat_c[i*matrix_sz+k] = CircuitExecution::circ_exec->xor_gate(mat_c[i*matrix_sz+k].bit, tmp.bit);
			}
		}
	}
	
    std::cout << HIGH64(mat_c[0].bit) << ' ' << LOW64(mat_c[0].bit) << std::endl;
    
	for (int i = 0; i < test_n; i++) {
        bool tmp = mat_c[i].reveal(PUBLIC);
        if (tmp == 1) exit(0);
    }
	bool cheated = finalize_zk_bool<BoolIO<NetIO>>();
	if(cheated) error("cheated\n");
	cout  << "\t" << time_from(start)<<" us\t"<<party<<endl;
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
