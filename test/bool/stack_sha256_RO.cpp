#include <emp-zk/emp-zk.h>
#include <iostream>
#include "emp-tool/emp-tool.h"
using namespace emp;
using namespace std;

int port, party;
const int threads = 1;
const string circuit_file_location = string("./sha256.txt");
const int and_gate_per_cir = 22573;

int gate_num;
int wire_num;
int input_size0, input_size1;
int output_size;
enum GType{
    ANDG,
    XORG,
    INVG
};
struct Bgate {
    int l, r;
    int o;
    GType gt;
};
vector<Bgate> bgate;
bool output_inv[256];

void read_circuit() {
    ifstream fin(circuit_file_location);
    fin >> gate_num >> wire_num;
    int n, m;
    fin >> n >> input_size0 >> input_size1;
    fin >> n >> output_size;
    bgate.resize(gate_num);    
    for (int i = 0; i < gate_num; i++) {
        fin >> n >> m;
        if (n == 2) fin >> bgate[i].l >> bgate[i].r;
        else fin >> bgate[i].l;
        fin >> bgate[i].o;
        string tmp;
        fin >> tmp;
        if (tmp[0] == 'A') bgate[i].gt = ANDG;
        else if (tmp[0] == 'X') bgate[i].gt = XORG;
        else if (tmp[0] == 'I') bgate[i].gt = INVG;
        else exit(-2);
    }
    fin.close();
}

void prepare_witness(int itr, bool *l, bool *r) {
    int andgateid = 0;
    bool wire[wire_num];
    for (int i = 0; i < input_size0 + input_size1; i++) wire[i] = 0;
    for (int it = 0; it < itr; it++) {
        for (int i = 0; i < gate_num; i++) {
            if (bgate[i].gt == ANDG) {
                wire[bgate[i].o] = wire[bgate[i].l] & wire[bgate[i].r];
                l[andgateid] = wire[bgate[i].l];
                r[andgateid++] = wire[bgate[i].r];
            } else if (bgate[i].gt == XORG) {
                wire[bgate[i].o] = wire[bgate[i].l] ^ wire[bgate[i].r];
            } else {
                wire[bgate[i].o] = !wire[bgate[i].l];
            }
        }
        for (int i = 0; i <  input_size1; i++) wire[input_size0 + i] = wire[wire_num - input_size1 + i];
    }
    for (int i = 0; i < output_size; i++) output_inv[i] = wire[wire_num-output_size+i];    
}

void mul_left_right(int itr, Bit *in, Bit *l, Bit *r, Bit *o, Bit &one, block *chis, block &authen, block &value) {

    int andgatesize = and_gate_per_cir * itr;

    int andgateid = 0;
    block wire[wire_num];
    for (int i = 0; i < input_size0 + input_size1; i++) wire[i] = in[i].bit;

    for (int it = 0; it < itr; it++) {
        for (int i = 0; i < gate_num; i++) {
            if (bgate[i].gt == ANDG) {
                //wire[bgate[i].o] = wire[bgate[i].l] & wire[bgate[i].r];
                //l[andgateid] = wire[bgate[i].l];
                //r[andgateid++] = wire[bgate[i].r];
                block tmp, tmp2;
                tmp = wire[bgate[i].l] ^ l[andgateid].bit;
                gfmul(tmp, chis[andgateid], &tmp2);
                authen = authen ^ tmp2;
                if (getLSB(tmp)) value = value ^ chis[andgateid];

                tmp = wire[bgate[i].r] ^ r[andgateid].bit;
                gfmul(tmp, chis[andgatesize + andgateid], &tmp2);
                authen = authen ^ tmp2;
                if (getLSB(tmp)) value = value ^ chis[andgatesize + andgateid];

                wire[bgate[i].o] = o[andgateid++].bit;
            } else if (bgate[i].gt == XORG) {
                wire[bgate[i].o] = wire[bgate[i].l] ^ wire[bgate[i].r];
            } else {
                wire[bgate[i].o] = wire[bgate[i].l] ^ one.bit;
            }
        }
        for (int i = 0; i <  input_size1; i++) wire[input_size0 + i] = wire[wire_num - input_size1 + i];
    }

    for (int i = 0; i < output_size; i++) {
        if (output_inv[i]) wire[input_size0 + i] = wire[input_size0 + i] ^ one.bit;

        block tmp2;
        gfmul(wire[input_size0 + i], chis[andgatesize*2 + i], &tmp2);
        authen = authen ^ tmp2;
        if (getLSB(wire[input_size0 + i])) value = value ^ chis[andgatesize * 2 + i];
    }

}


void test_circuit_zk(BoolIO<NetIO> *ios[threads], int party, int rep_sz, int branch_sz) {

    read_circuit();

    bool *lvalue = new bool[and_gate_per_cir * rep_sz];
    bool *rvalue = new bool[and_gate_per_cir * rep_sz];
    auto st = clock_start();
    prepare_witness(rep_sz, lvalue, rvalue);
    cout  << "\t" << time_from(st)<<" us\t"<<party<<endl;

    
    setup_zk_bool<BoolIO<NetIO>>(ios, threads, party);
    Bit *in = new Bit[input_size0+input_size1];
    Bit *l = new Bit[and_gate_per_cir * rep_sz];
    Bit *r = new Bit[and_gate_per_cir * rep_sz];
    Bit *o = new Bit[and_gate_per_cir * rep_sz];
    int chisize = and_gate_per_cir * rep_sz * 2 + output_size;
    block *chis = new block[chisize];    

    cout  << "\t" << time_from(st)<<" us\t"<<party<<endl;
    cout << "HAHHA" << endl << endl;

	
	auto start = clock_start();

    for (int i = 0; i < input_size0+input_size1; i++) in[i] = Bit(false, ALICE);    
    for (int i = 0; i < and_gate_per_cir * rep_sz; i++) {
        l[i] = Bit(lvalue[i], ALICE);
        r[i] = Bit(rvalue[i], ALICE);
        o[i] = l[i] & r[i];
    }
    //delete[] lvalue;
    //delete[] rvalue;    

    block s_seed; 
    if (party == ALICE) {
		ios[0]->recv_data(&s_seed, sizeof(block));
    } else {
        PRG().random_block(&s_seed, 1);
        ios[0]->send_data(&s_seed, sizeof(block));
        ios[0]->flush();
    }
    PRG prg_s(&s_seed);       

    Bit one = Bit(true, PUBLIC);

	//block chi; // random challenge
	block proofs[branch_sz];
    block values[branch_sz];
    for (int i = 0; i < branch_sz; i++) {
        proofs[i] = makeBlock(0, 0);
        values[i] = makeBlock(0, 0);
    }

	if (party == ALICE) {		
        /*
        uint64_t tmp1, tmp2;
        ios[0]->recv_data(&tmp1, sizeof(uint64_t));
        ios[0]->recv_data(&tmp2, sizeof(uint64_t));
        chi = makeBlock(tmp1, tmp2);

        chis[0] = makeBlock(0, 1);*/
        for (int i = 0; i < chisize; i++) prg_s.random_block(&chis[i], 1);
        //gfmul(chis[i-1], chi, &chis[i]);

        for (int i = 0; i < branch_sz; i++) mul_left_right(rep_sz, in, l, r, o, one, chis, proofs[i], values[i]);
    } else {
        /*
		PRG tmpprg;
        uint64_t tmp1, tmp2;
		tmpprg.random_data(&tmp1, sizeof(uint64_t));
        tmpprg.random_data(&tmp2, sizeof(uint64_t));
        ios[0]->send_data(&tmp1, sizeof(uint64_t));
        ios[0]->send_data(&tmp2, sizeof(uint64_t));
        ios[0]->flush();
        chi = makeBlock(tmp1, tmp2);

        chis[0] = makeBlock(0, 1);*/
        for (int i = 0; i < chisize; i++) prg_s.random_block(&chis[i], 1);
        //gfmul(chis[i-1], chi, &chis[i]);        

        for (int i = 0; i < branch_sz; i++) mul_left_right(rep_sz, in, l, r, o, one, chis, proofs[i], values[i]);        
    }
    //delete[] chis;

    // for (int i = 0; i < branch_sz; i++) cout << HIGH64(proofs[i]) << ' ' << LOW64(proofs[i]) << ' ' << HIGH64(values[i]) << ' ' << LOW64(values[i]) << std::endl;
    // final proof

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

	if(argc < 6) {
		std::cout << "usage: a.out PARTY(1/2) PORT ADDR #ITERATION #BRANCH" << std::endl;
		return -1;
	}	
	
	parse_party_and_port(argv, &party, &port);
	BoolIO<NetIO>* ios[threads];
	for(int i = 0; i < threads; ++i)
		ios[i] = new BoolIO<NetIO>(new NetIO(party == ALICE?nullptr:argv[3],port+i), party==ALICE);

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
		num = atoi(argv[4]);
        branch = atoi(argv[5]);
	}

	test_circuit_zk(ios, party, num, branch);

	for(int i = 0; i < threads; ++i) {
		delete ios[i]->io;
		delete ios[i];
	}
	return 0;
}
