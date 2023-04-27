#include <emp-zk/emp-zk.h>
#include <iostream>
#include "emp-tool/emp-tool.h"
using namespace emp;
using namespace std;

int port, party;
const int threads = 1;

void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int matrix_sz, int branch_sz) {



	long long test_n = matrix_sz * matrix_sz;
	long long mul_sz = matrix_sz * matrix_sz * matrix_sz;

	setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);

    block *chis = new block[mul_sz*2+test_n];
	Bit *mat_a = new Bit[test_n];
	Bit *mat_b = new Bit[test_n];
	Bit *mul_le = new Bit[mul_sz];
	Bit *mul_ri = new Bit[mul_sz];
	Bit *mul_ou = new Bit[mul_sz];

	auto start = clock_start();


	// Input
	for (int i = 0; i < test_n; i++) {
		mat_a[i] = Bit(true, ALICE);
		mat_b[i] = Bit(true, ALICE);
	}

	// Left, Right, Output of AND gates
	for (int i = 0; i < mul_sz; i++) {
		mul_le[i] = Bit(true, ALICE);
		mul_ri[i] = Bit(true, ALICE);
		mul_ou[i] = mul_le[i] & mul_ri[i];
	}

	// Unit for constant offsets
	// IntFp one = IntFp(1, PUBLIC);
    Bit one = Bit(true, PUBLIC);
	block chi; // random challenge
	//uint64_t delta = ZKFpExec::zk_exec->get_delta();
	//std::cout << "DELTA = " << delta << std::endl;

    //block *chis = new block[];

	block proofs[branch_sz];
    block values[branch_sz];

	if (party == ALICE) {		
        uint64_t tmp1, tmp2;
        ios[0]->recv_data(&tmp1, sizeof(uint64_t));
        ios[0]->recv_data(&tmp2, sizeof(uint64_t));
        chi = makeBlock(tmp1, tmp2);

        //block *chis = new block[mul_sz*2+test_n];
        chis[0] = makeBlock(0, 1);
        for (int i = 1; i < mul_sz*2 + test_n; i++) gfmul(chis[i-1], chi, &chis[i]);

		// Bob needs to generate potential proofs for all branches
		for (int bid = 0; bid < branch_sz; bid++) {

            proofs[bid] = makeBlock(0, 0);
            values[bid] = makeBlock(0, 0);
			//block coeff = makeBlock(0, 1); // coeff as chi^k			
            uint64_t coe_id = 0;

			// First 2*|C_MUL| rows to check that le/ri are formed correctly
			// Checking left
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id / matrix_sz;

                block tmpmul;
                gfmul(chis[coe_id], mul_le[mul_id].bit ^ mat_a[i].bit, &tmpmul);
                proofs[bid] = proofs[bid] ^ tmpmul;

                bool tmp = getLSB(mul_le[mul_id].bit) ^ getLSB(mat_a[i].bit);
                if (tmp) values[bid] = values[bid] ^ chis[coe_id];

                //block newcoeff;
                //gfmul(coeff, chi, &newcoeff);
                //coeff = newcoeff;
                coe_id++;
			}
			// Checking right
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id % test_n;

                block tmpmul;
                gfmul(chis[coe_id], mul_ri[mul_id].bit ^ mat_b[i].bit, &tmpmul);
                proofs[bid] = proofs[bid] ^ tmpmul;

                bool tmp = getLSB(mul_ri[mul_id].bit) ^ getLSB(mat_b[i].bit);                
                if (tmp) values[bid] = values[bid] ^ chis[coe_id];

                //block newcoeff;
                //gfmul(coeff, chi, &newcoeff);
                //coeff = newcoeff;
                coe_id++;
			}	
			// Checking output: the active branch is the one that res[i][j] = matrix_sz
			for (int i = 0; i < matrix_sz; i++) 
				for (int k = 0; k < matrix_sz; k++) {

                    block aut_tmp = one.bit;
                    bool val_tmp = true;

                    if (bid % 2 == 0) {
                        aut_tmp = makeBlock(0, 0);
                        val_tmp = false;
                    }

					for (int j = 0; j < matrix_sz; j++) {
                        aut_tmp = aut_tmp ^ mul_ou[i*test_n + j*matrix_sz + k].bit;
                        val_tmp = val_tmp ^ getLSB(mul_ou[i*test_n + j*matrix_sz + k].bit);
					}

                    block tmpmul;
                    gfmul(chis[coe_id], aut_tmp, &tmpmul);
                    proofs[bid] = proofs[bid] ^ tmpmul;

                    if (val_tmp) values[bid] = values[bid] ^ chis[coe_id];

                    //block newcoeff;
                    //gfmul(coeff, chi, &newcoeff);
                    //coeff = newcoeff;			
                    coe_id++;
				}							
		}

        //delete[] chis;		

	} else {

		PRG tmpprg;
        uint64_t tmp1, tmp2;
		tmpprg.random_data(&tmp1, sizeof(uint64_t));
        tmpprg.random_data(&tmp2, sizeof(uint64_t));
        ios[0]->send_data(&tmp1, sizeof(uint64_t));
        ios[0]->send_data(&tmp2, sizeof(uint64_t));
        ios[0]->flush();
        chi = makeBlock(tmp1, tmp2);

        
        chis[0] = makeBlock(0, 1);
        for (int i = 1; i < mul_sz*2 + test_n; i++) gfmul(chis[i-1], chi, &chis[i]);
		
		// Bob needs to generate potential proofs for all branches
		for (int bid = 0; bid < branch_sz; bid++) {

            proofs[bid] = makeBlock(0, 0);
			// block coeff = makeBlock(0, 1); // coeff as chi^k			
            uint64_t coe_id = 0;

			// First 2*|C_MUL| rows to check that le/ri are formed correctly
			// Checking left
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id / matrix_sz;
                block tmpmul;
                gfmul(chis[coe_id], mul_le[mul_id].bit ^ mat_a[i].bit, &tmpmul);
                proofs[bid] = proofs[bid] ^ tmpmul;

                //block newcoeff;
                //gfmul(coeff, chi, &newcoeff);
                //coeff = newcoeff;
                coe_id++;
			}
			// Checking right
			for (int mul_id = 0; mul_id < mul_sz; mul_id++) {
				int i = mul_id % test_n;

                block tmpmul;
                gfmul(chis[coe_id], mul_ri[mul_id].bit ^ mat_b[i].bit, &tmpmul);
                proofs[bid] = proofs[bid] ^ tmpmul;

                //block newcoeff;
                //gfmul(coeff, chi, &newcoeff);
                //coeff = newcoeff;
                coe_id++;
			}	
			// Checking output: the active branch is the one that res[i][j] = matrix_sz
			for (int i = 0; i < matrix_sz; i++) 
				for (int k = 0; k < matrix_sz; k++) {

                    block aut_tmp = one.bit;

                    if (bid % 2 == 0) {
                        aut_tmp = makeBlock(0, 0);
                    }

					for (int j = 0; j < matrix_sz; j++) {
                        aut_tmp = aut_tmp ^ mul_ou[i*test_n + j*matrix_sz + k].bit;
					}

                    block tmpmul;
                    gfmul(chis[coe_id], aut_tmp, &tmpmul);
                    proofs[bid] = proofs[bid] ^ tmpmul;

                    //block newcoeff;
                    //gfmul(coeff, chi, &newcoeff);
                    //coeff = newcoeff;			
                    coe_id++;
				}							
		}
        //delete[] chis;
	}	

    // final proof
    //for (int i = 0; i < branch_sz; i++) cout << HIGH64(proofs[i]) << LOW64(proofs[i]) << endl;
    //for (int i = 0; i < branch_sz; i++) cout << HIGH64(values[i]) << LOW64(values[i]) << endl;

    block aut[branch_sz], val[branch_sz];
    for (int i = 0; i < branch_sz; i++)
        zkp_get_ope<BoolIO<NetIO>>(aut[i], val[i]);
    block delta = get_bool_delta<BoolIO<NetIO>>(party);

    if (party == ALICE) {
        // commit the intermiate wires of final MUL gates
        block inter_mul[branch_sz];
        gfmul(values[0], values[1], &inter_mul[0]);
        for (int i = 1; i < branch_sz-1; i++) gfmul(inter_mul[i-1], values[i+1], &inter_mul[i]);
        for (int i = 0; i < branch_sz-1; i++) {
            val[i] = val[i] ^ inter_mul[i];
            ios[0]->send_data(&val[i], sizeof(block));
            val[i] = inter_mul[i];
        }
        ios[0]->flush();

        // LPZK
        uint64_t tmp1, tmp2;
        ios[0]->recv_data(&tmp1, sizeof(uint64_t));
        ios[0]->recv_data(&tmp2, sizeof(uint64_t));
        block kai = makeBlock(tmp1, tmp2);

        block coeff = makeBlock(0, 1);
        block A0 = makeBlock(0, 0);
        block A1 = makeBlock(0, 0);

        // first mul gate
        block tmp;
        gfmul(proofs[0], proofs[1] ,&tmp);
        block tmpmul;
        gfmul(coeff, tmp, &tmpmul);
        A0 = A0 ^ tmpmul;

        gfmul(values[0], proofs[1], &tmp);
        gfmul(coeff, tmp, &tmpmul);
        A1 = A1 ^ tmpmul;
        gfmul(values[1], proofs[0], &tmp);
        gfmul(coeff, tmp, &tmpmul);
        A1 = A1 ^ tmpmul;
        gfmul(coeff, aut[0], &tmpmul);
        A1 = A1 ^ tmpmul;

        gfmul(kai, coeff, &tmpmul);
        coeff = tmpmul;

        for (int i = 1; i < branch_sz-1; i++) {            
            gfmul(aut[i-1], proofs[i+1], &tmp);
            gfmul(coeff, tmp, &tmpmul);
            A0 = A0 ^ tmpmul;

            gfmul(val[i-1], proofs[i+1], &tmp);
            gfmul(coeff, tmp, &tmpmul);
            A1 = A1 ^ tmpmul;
            gfmul(aut[i-1], values[i+1], &tmp);            
            gfmul(coeff, tmp, &tmpmul);
            A1 = A1 ^ tmpmul;
            gfmul(coeff, aut[i], &tmpmul);
            A1 = A1 ^ tmpmul;

            gfmul(kai, coeff, &tmpmul);
            coeff = tmpmul;
        }

        // add OTP for LPZK
        A0 = A0 ^ aut[branch_sz-1];
        A1 = A1 ^ val[branch_sz-1];
        ios[0]->send_data(&A0, sizeof(block));
        ios[0]->send_data(&A1, sizeof(block));
        //final zero        
        ios[0]->send_data(&aut[branch_sz-2], sizeof(block));
        ios[0]->flush();

    } else {
        for (int i = 0; i < branch_sz-1; i++) {
            block tmp;
            ios[0]->recv_data(&tmp, sizeof(block));
            block tmpmul;
            gfmul(delta, tmp, &tmpmul);
            aut[i] = aut[i] ^ tmpmul; 
        }

        // LPZK
		PRG tmpprg;
        uint64_t tmp1, tmp2;
		tmpprg.random_data(&tmp1, sizeof(uint64_t));
        tmpprg.random_data(&tmp2, sizeof(uint64_t));
        ios[0]->send_data(&tmp1, sizeof(uint64_t));
        ios[0]->send_data(&tmp2, sizeof(uint64_t));
        ios[0]->flush();
        block kai = makeBlock(tmp1, tmp2);

        block coeff = makeBlock(0, 1);
        block accum = makeBlock(0, 0);

        // first mul gate
        block tmp;
        gfmul(proofs[0], proofs[1] ,&tmp);
        block tmpmul;
        gfmul(coeff, tmp, &tmpmul);
        accum = accum ^ tmpmul;
        gfmul(aut[0], delta, &tmp);
        gfmul(coeff, tmp, &tmpmul);
        accum = accum ^ tmpmul;
        gfmul(kai, coeff, &tmpmul);
        coeff = tmpmul;

        for (int i = 1; i < branch_sz-1; i++) {
            gfmul(aut[i-1], proofs[i+1], &tmp);
            gfmul(coeff, tmp, &tmpmul);
            accum = accum ^ tmpmul;
            gfmul(aut[i], delta, &tmp);
            gfmul(coeff, tmp, &tmpmul);
            accum = accum ^ tmpmul;
            gfmul(kai, coeff, &tmpmul);
            coeff = tmpmul;
        }
        accum = accum ^ aut[branch_sz-1];

        block A0, A1;
        ios[0]->recv_data(&A0, sizeof(block));
        ios[0]->recv_data(&A1, sizeof(block));
        gfmul(A1, delta, &tmp);
        tmp = tmp ^ A0;

        if (HIGH64(tmp) == HIGH64(accum) && LOW64(tmp) == LOW64(accum)) cout << "FINAL LPZK SUCCESS" << endl;
        else error("cheated\n");

        // final 0 check
        ios[0]->recv_data(&tmp, sizeof(block));
        if (HIGH64(tmp) == HIGH64(aut[branch_sz-2]) && LOW64(tmp) == LOW64(aut[branch_sz-2])) cout << "FINAL 0 SUCCESS" << endl;
        else error("cheated\n");

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
