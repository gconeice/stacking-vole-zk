#ifndef GFMUL_H__
#define GFMUL_H__

#include "emp-tool/emp-tool.h"
#include <iostream>

template<int N>
inline void mul128_0(__m128i *a, __m128i *b, __m128i *res1, __m128i *res2) {
	__m128i tmp3[N], tmp4[N], tmp5[N], tmp6[N];
	for(int i = 0; i < N; ++i) { 
		tmp3[i] = _mm_clmulepi64_si128(a[i], b[i], 0x00);
		tmp4[i] = _mm_clmulepi64_si128(a[i], b[i], 0x10);
		tmp5[i] = _mm_clmulepi64_si128(a[i], b[i], 0x01);
		tmp6[i] = _mm_clmulepi64_si128(a[i], b[i], 0x11);
	}
	for(int i = 0; i < N; ++i) {
		tmp4[i] = _mm_xor_si128(tmp4[i], tmp5[i]);
		tmp5[i] = _mm_slli_si128(tmp4[i], 8);
		tmp4[i] = _mm_srli_si128(tmp4[i], 8);
		tmp3[i] = _mm_xor_si128(tmp3[i], tmp5[i]);
		tmp6[i] = _mm_xor_si128(tmp6[i], tmp4[i]);
		res1[i] = tmp3[i];
		res2[i] = tmp6[i];
	}
}

inline void reduction0(__m128i tmp3, __m128i tmp6, __m128i* res) {
	/*__m128i tmp7, tmp8, tmp9, tmp10, tmp11, tmp12;
	__m128i XMMMASK = _mm_setr_epi32(0xffffffff, 0x0, 0x0, 0x0);
	tmp7 = _mm_srli_epi32(tmp6, 31);
	tmp8 = _mm_srli_epi32(tmp6, 30);
	tmp9 = _mm_srli_epi32(tmp6, 25);
	tmp7 = _mm_xor_si128(tmp7, tmp8);
	tmp7 = _mm_xor_si128(tmp7, tmp9);
	tmp8 = _mm_shuffle_epi32(tmp7, 147);

	tmp7 = _mm_and_si128(XMMMASK, tmp8);
	tmp8 = _mm_andnot_si128(XMMMASK, tmp8);
	tmp3 = _mm_xor_si128(tmp3, tmp8);
	tmp6 = _mm_xor_si128(tmp6, tmp7);
	tmp10 = _mm_slli_epi32(tmp6, 1);
	tmp3 = _mm_xor_si128(tmp3, tmp10);
	tmp11 = _mm_slli_epi32(tmp6, 2);
	tmp3 = _mm_xor_si128(tmp3, tmp11);
	tmp12 = _mm_slli_epi32(tmp6, 7);
	tmp3 = _mm_xor_si128(tmp3, tmp12);

	*res = _mm_xor_si128(tmp3, tmp6);*/
    __m128i tmp2, tmp4, tmp5, tmp7, tmp8, tmp9;
    tmp7 = _mm_srli_epi32(tmp3, 31);
    tmp8 = _mm_srli_epi32(tmp6, 31);
    tmp3 = _mm_slli_epi32(tmp3, 1);
    tmp6 = _mm_slli_epi32(tmp6, 1);

    tmp9 = _mm_srli_si128(tmp7, 12);
    tmp8 = _mm_slli_si128(tmp8, 4);
    tmp7 = _mm_slli_si128(tmp7, 4);
    tmp3 = _mm_or_si128(tmp3, tmp7);
    tmp6 = _mm_or_si128(tmp6, tmp8);
    tmp6 = _mm_or_si128(tmp6, tmp9);

    tmp7 = _mm_slli_epi32(tmp3, 31);
    tmp8 = _mm_slli_epi32(tmp3, 30);
    tmp9 = _mm_slli_epi32(tmp3, 25);
    tmp7 = _mm_xor_si128(tmp7, tmp8);
    tmp7 = _mm_xor_si128(tmp7, tmp9);
    tmp8 = _mm_srli_si128(tmp7, 4);
    tmp7 = _mm_slli_si128(tmp7, 12);
    tmp3 = _mm_xor_si128(tmp3, tmp7);

    tmp2 = _mm_srli_epi32(tmp3, 1);
    tmp4 = _mm_srli_epi32(tmp3, 2);
    tmp5 = _mm_srli_epi32(tmp3, 7);
    tmp2 = _mm_xor_si128(tmp2, tmp4);
    tmp2 = _mm_xor_si128(tmp2, tmp5);
    tmp2 = _mm_xor_si128(tmp2, tmp8);
    tmp3 = _mm_xor_si128(tmp3, tmp2);
    *res = _mm_xor_si128(tmp6, tmp3);
}

template<int N>
inline void gfmul_1(__m128i *a, __m128i *b, __m128i *res) {
	__m128i tmp3[N], tmp6[N];
	mul128_0<N>(a, b, tmp3, tmp6);
	for(int i = 0; i < N; ++i)
		reduction0(tmp3[i], tmp6[i], res+i);
}

/* batch 3 */
inline void mul128_3(__m128i &a1, __m128i &a2, __m128i &a3,
		__m128i &b1, __m128i &b2, __m128i &b3,
		__m128i *res1, __m128i *res2) {
	__m128i tmp3[3], tmp4[3], tmp5[3], tmp6[3];
	tmp3[0] = _mm_clmulepi64_si128(a1, b1, 0x00);
	tmp4[0] = _mm_clmulepi64_si128(a1, b1, 0x10);
	tmp5[0] = _mm_clmulepi64_si128(a1, b1, 0x01);
	tmp6[0] = _mm_clmulepi64_si128(a1, b1, 0x11);
	tmp3[1] = _mm_clmulepi64_si128(a2, b2, 0x00);
	tmp4[1] = _mm_clmulepi64_si128(a2, b2, 0x10);
	tmp5[1] = _mm_clmulepi64_si128(a2, b2, 0x01);
	tmp6[1] = _mm_clmulepi64_si128(a2, b2, 0x11);
	tmp3[2] = _mm_clmulepi64_si128(a3, b3, 0x00);
	tmp4[2] = _mm_clmulepi64_si128(a3, b3, 0x10);
	tmp5[2] = _mm_clmulepi64_si128(a3, b3, 0x01);
	tmp6[2] = _mm_clmulepi64_si128(a3, b3, 0x11);
	for(int i = 0; i < 3; ++i) {
		tmp4[i] = _mm_xor_si128(tmp4[i], tmp5[i]);
		tmp5[i] = _mm_slli_si128(tmp4[i], 8);
		tmp4[i] = _mm_srli_si128(tmp4[i], 8);
		tmp3[i] = _mm_xor_si128(tmp3[i], tmp5[i]);
		tmp6[i] = _mm_xor_si128(tmp6[i], tmp4[i]);
		res1[i] = tmp3[i];
		res2[i] = tmp6[i];
	}
}

inline void gfmul3(__m128i &a1, __m128i &a2, __m128i &a3,
		__m128i &b1, __m128i &b2, __m128i &b3,
	       	__m128i *res) {
	__m128i tmp3[3], tmp6[3];
	mul128_3(a1, a2, a3, b1, b2, b3, tmp3, tmp6);
	for(int i = 0; i < 3; ++i)
		reduction0(tmp3[i], tmp6[i], res+i);
}


/* batch 4 */
inline void mul128_4(__m128i &a1, __m128i &a2, __m128i &a3, __m128i &a4,
		__m128i &b1, __m128i &b2, __m128i &b3, __m128i &b4,
		__m128i *res1, __m128i *res2) {
	__m128i tmp3[4], tmp4[4], tmp5[4], tmp6[4];
	tmp3[0] = _mm_clmulepi64_si128(a1, b1, 0x00);
	tmp4[0] = _mm_clmulepi64_si128(a1, b1, 0x10);
	tmp5[0] = _mm_clmulepi64_si128(a1, b1, 0x01);
	tmp6[0] = _mm_clmulepi64_si128(a1, b1, 0x11);
	tmp3[1] = _mm_clmulepi64_si128(a2, b2, 0x00);
	tmp4[1] = _mm_clmulepi64_si128(a2, b2, 0x10);
	tmp5[1] = _mm_clmulepi64_si128(a2, b2, 0x01);
	tmp6[1] = _mm_clmulepi64_si128(a2, b2, 0x11);
	tmp3[2] = _mm_clmulepi64_si128(a3, b3, 0x00);
	tmp4[2] = _mm_clmulepi64_si128(a3, b3, 0x10);
	tmp5[2] = _mm_clmulepi64_si128(a3, b3, 0x01);
	tmp6[2] = _mm_clmulepi64_si128(a3, b3, 0x11);
	tmp3[3] = _mm_clmulepi64_si128(a4, b4, 0x00);
	tmp4[3] = _mm_clmulepi64_si128(a4, b4, 0x10);
	tmp5[3] = _mm_clmulepi64_si128(a4, b4, 0x01);
	tmp6[3] = _mm_clmulepi64_si128(a4, b4, 0x11);
	for(int i = 0; i < 4; ++i) {
		tmp4[i] = _mm_xor_si128(tmp4[i], tmp5[i]);
		tmp5[i] = _mm_slli_si128(tmp4[i], 8);
		tmp4[i] = _mm_srli_si128(tmp4[i], 8);
		tmp3[i] = _mm_xor_si128(tmp3[i], tmp5[i]);
		tmp6[i] = _mm_xor_si128(tmp6[i], tmp4[i]);
		res1[i] = tmp3[i];
		res2[i] = tmp6[i];
	}
}

inline void gfmul4(__m128i &a1, __m128i &a2, __m128i &a3, __m128i &a4,
		__m128i &b1, __m128i &b2, __m128i &b3, __m128i &b4,
	       	__m128i *res) {
	__m128i tmp3[4], tmp6[4];
	mul128_4(a1, a2, a3, a4, b1, b2, b3, b4, tmp3, tmp6);
	for(int i = 0; i < 4; ++i)
		reduction0(tmp3[i], tmp6[i], res+i);
}


/* batch 6 */
inline void mul128_6(__m128i &a1, __m128i &a2, __m128i &a3,
		__m128i &a4, __m128i &a5, __m128i &a6,
		__m128i &b1, __m128i &b2, __m128i &b3,
		__m128i &b4, __m128i &b5, __m128i &b6,
		__m128i *res1, __m128i *res2) {
	__m128i tmp3[6], tmp4[6], tmp5[6], tmp6[6];
	tmp3[0] = _mm_clmulepi64_si128(a1, b1, 0x00);
	tmp4[0] = _mm_clmulepi64_si128(a1, b1, 0x10);
	tmp5[0] = _mm_clmulepi64_si128(a1, b1, 0x01);
	tmp6[0] = _mm_clmulepi64_si128(a1, b1, 0x11);
	tmp3[1] = _mm_clmulepi64_si128(a2, b2, 0x00);
	tmp4[1] = _mm_clmulepi64_si128(a2, b2, 0x10);
	tmp5[1] = _mm_clmulepi64_si128(a2, b2, 0x01);
	tmp6[1] = _mm_clmulepi64_si128(a2, b2, 0x11);
	tmp3[2] = _mm_clmulepi64_si128(a3, b3, 0x00);
	tmp4[2] = _mm_clmulepi64_si128(a3, b3, 0x10);
	tmp5[2] = _mm_clmulepi64_si128(a3, b3, 0x01);
	tmp6[2] = _mm_clmulepi64_si128(a3, b3, 0x11);
	tmp3[3] = _mm_clmulepi64_si128(a4, b4, 0x00);
	tmp4[3] = _mm_clmulepi64_si128(a4, b4, 0x10);
	tmp5[3] = _mm_clmulepi64_si128(a4, b4, 0x01);
	tmp6[3] = _mm_clmulepi64_si128(a4, b4, 0x11);
	tmp3[4] = _mm_clmulepi64_si128(a5, b5, 0x00);
	tmp4[4] = _mm_clmulepi64_si128(a5, b5, 0x10);
	tmp5[4] = _mm_clmulepi64_si128(a5, b5, 0x01);
	tmp6[4] = _mm_clmulepi64_si128(a5, b5, 0x11);
	tmp3[5] = _mm_clmulepi64_si128(a6, b6, 0x00);
	tmp4[5] = _mm_clmulepi64_si128(a6, b6, 0x10);
	tmp5[5] = _mm_clmulepi64_si128(a6, b6, 0x01);
	tmp6[5] = _mm_clmulepi64_si128(a6, b6, 0x11);
	for(int i = 0; i < 6; ++i) {
		tmp4[i] = _mm_xor_si128(tmp4[i], tmp5[i]);
		tmp5[i] = _mm_slli_si128(tmp4[i], 8);
		tmp4[i] = _mm_srli_si128(tmp4[i], 8);
		tmp3[i] = _mm_xor_si128(tmp3[i], tmp5[i]);
		tmp6[i] = _mm_xor_si128(tmp6[i], tmp4[i]);
		res1[i] = tmp3[i];
		res2[i] = tmp6[i];
	}
}

inline void gfmul6(__m128i &a1, __m128i &a2, __m128i &a3,
		__m128i &a4, __m128i &a5, __m128i &a6,
		__m128i &b1, __m128i &b2, __m128i &b3,
		__m128i &b4, __m128i &b5, __m128i &b6,
		__m128i *res) {
	__m128i tmp3[6], tmp6[6];
	mul128_6(a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6, tmp3, tmp6);
	for(int i = 0; i < 6; ++i)
		reduction0(tmp3[i], tmp6[i], res+i);
}
#endif

/*int main() {
	PRG prg;
	block res[8];
	block res2[8];
	block res3[8];
	prg.random_block(res, 4);
	prg.random_block(res2, 4);
	prg.random_block(res3, 4);

	cout << "Benchmark: ";
	auto start = clock_start();
	for(int i = 0; i < 1024*1024*32; ++i) {
		gfmul_1(res, res2, res3);
		gfmul_1(res3, res, res);
	}
	cout << 1024*1024*32*2/(time_from(start))*1e6 << " operations/second" << endl;
	cout <<res[0]<<endl;
	return 0;
}*/
