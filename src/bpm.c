
#include "bpm.h"
#include  <stdalign.h>

#ifdef HAVE_AVX2
#include <immintrin.h>

#define MACRO_MIN(a,b)          (((a)<(b))?(a):(b))

__m256i BROADCAST_MASK[16];

 void bitShiftLeft256ymm (__m256i *data, int count);
__m256i bitShiftRight256ymm (__m256i *data, int count);

/* taken from Alexander Yee: http://www.numberworld.org/y-cruncher/internals/addition.html#ks_add */
 __m256i add256(uint32_t carry, __m256i A, __m256i B);
#endif


uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m)
{
        register uint64_t VP,VN,D0,HN,HP,X;
        register uint64_t i;//,c;
        uint64_t MASK = 0;
        int64_t diff;
        uint64_t B[13];
        int8_t k;

        if(m > 63){
                m = 63;
        }
        diff = m;
        k = m;
        for(i = 0; i < 13;i++){
                B[i] = 0;
        }

        for(i = 0; i < m;i++){
                B[p[i]] |= (1ul << i);
        }

        VP = (1ul << (m))-1 ;

        VN = 0ul;

        m--;
        MASK = 1ul << (m);

        //fprintf(stdout,"BEGINNING\t%lu %lu\n",VN,VP);
        for(i = 0; i < n;i++){
                //        fprintf(stdout,"%lu:\t",i);
                X = (B[t[i]] | VN);
                //fprintf(stdout,"%lu ", X);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                //fprintf(stdout,"%lu ", D0);
                HN = VP & D0;
                //fprintf(stdout,"%lu ", HN);
                HP = VN | (~(VP | D0));
                //fprintf(stdout,"%lu ", HP);
                X = HP << 1ul;
                //fprintf(stdout,"%lu ", X);
                VN = X & D0;
                //fprintf(stdout,"%lu ", VN);

                VP = (HN << 1ul) | (~(X | D0));
                //fprintf(stdout,"%lu ", VP);


                diff += (HP & MASK)? 1 : 0;// >> m;
                diff -= (HN & MASK)? 1 : 0;// >> m;
                if(diff < k){
                        k = diff;

                }
                //fprintf(stdout,"%lu ", diff);
                //fprintf(stdout,"\n");

        }
        return k;
}


#ifdef HAVE_AVX2
uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m)
{
        __m256i VP,VN,D0,HN,HP,X,NOTONE;
        __m256i xmm1,xmm2;
        __m256i MASK;
        __m256i B[13];

        int i,j, k,diff;

        alignas(32)  uint32_t f[21][8];
        //int ALIGNED_(64) f[8];
        if(m > 255){
                m = 255;
        }

        for(i = 0; i < 21;i++){
                for(j = 0;j < 8;j++){
                        f[i][j] =0u;
                }
        }

        for(i = 0; i < m;i++){
                f[p[i]][i/32] |= (1 << (i % 32));
        }

        for(i = 0; i < 13;i++){
                B[i] = _mm256_load_si256((__m256i const*) &f[i]);
        }

        diff = m;
        k = m;

        VP     = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFul);
        VN     = _mm256_setzero_si256();
        NOTONE = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFFul);
        MASK   = _mm256_set_epi64x (0ul,0ul,0ul,1);
        m--;

        i = m / 64;
        while(i){
                bitShiftLeft256ymm(&MASK,64);
                i--;
        }
        bitShiftLeft256ymm(&MASK,m%64);

        for(i = 0; i < n ;i++){
                //X = (B[(int) *t] | VN);

                X = _mm256_or_si256(B[t[i]], VN);
                                //fprintf(stdout,"%lu ", X);
                //D0 = ((VP+(X&VP)) ^ VP) | X ;

                //print_256(X);
                xmm1 = _mm256_and_si256(X, VP);

                xmm2 = add256(0, VP, xmm1);
                //xmm2 = _mm256_add_epi64(VP, xmm1);
                xmm1 = _mm256_xor_si256(xmm2, VP);
                D0 = _mm256_or_si256(xmm1, X);
                //print_256(D0);
                //HN = VP & D0;
                HN =_mm256_and_si256(VP, D0);
                //print_256(HN);
                //HP = VN | ~(VP | D0);

                xmm1 = _mm256_or_si256(VP, D0);
                xmm2 = _mm256_andnot_si256(xmm1, NOTONE);
                HP = _mm256_or_si256(VN, xmm2);
                //print_256(HP);

                //X = HP << 1ul;
                X = HP;
                bitShiftLeft256ymm(&X, 1);

                //print_256(X);
                //VN = X & D0;
                VN= _mm256_and_si256(X, D0);
                //print_256(VN);
                //VP = (HN << 1ul) | ~(X | D0);
                xmm1 = HN;
                bitShiftLeft256ymm(&xmm1, 1);

                xmm2 = _mm256_or_si256(X, D0);

                xmm2 = _mm256_andnot_si256(xmm2, NOTONE);
                //xmm2 = _mm_andnot_si128 (xmm2,NOTONE);
                VP = _mm256_or_si256(xmm1, xmm2);
                //print_256(VP);


                //diff += (HP & MASK) >> m;
                diff += 1- _mm256_testz_si256(HP, MASK);

                ///diff -= (HN & MASK) >> m;
                diff -= 1- _mm256_testz_si256(HN,MASK);

                //fprintf(stdout,"%d ",diff);
                //xmm1 = _mm256_cmpgt_epi64(K, diff);
                k = MACRO_MIN(k, diff);
        }
        return k;
}




/* Must be called before BPM_256 is!!!  */
void set_broadcast_mask(void)
{
        BROADCAST_MASK[0] =  _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[1] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[2] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[3] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[4] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[5] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[6] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[7] = _mm256_set_epi64x(0x8000000000000000, 0x8000000000000001, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[8] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[9] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[10] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[11] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000000, 0x8000000000000001, 0x8000000000000001);
        BROADCAST_MASK[12] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000000, 0x8000000000000000);
        BROADCAST_MASK[13] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000000, 0x8000000000000001);
        BROADCAST_MASK[14] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000001, 0x8000000000000000);
        BROADCAST_MASK[15] = _mm256_set_epi64x(0x8000000000000001, 0x8000000000000001, 0x8000000000000001, 0x8000000000000001);
}


__m256i add256(uint32_t carry, __m256i A, __m256i B)
{
        A = _mm256_xor_si256(A, _mm256_set1_epi64x(0x8000000000000000));
        __m256i s = _mm256_add_epi64(A, B);
        __m256i cv = _mm256_cmpgt_epi64(A, s);
        __m256i mv = _mm256_cmpeq_epi64(s, _mm256_set1_epi64x(0x7fffffffffffffff));
        uint32_t c = _mm256_movemask_pd(_mm256_castsi256_pd(cv));
        uint32_t m = _mm256_movemask_pd(_mm256_castsi256_pd(mv));

        {
                c = m + 2*c; //  lea
                carry += c;
                m ^= carry;
                carry >>= 4;
                m &= 0x0f;
        }
        return _mm256_add_epi64(s, BROADCAST_MASK[m]);
}


//----------------------------------------------------------------------------
// bit shift left a 256-bit value using ymm registers
//          __m256i *data - data to shift
//          int count     - number of bits to shift
// return:  __m256i       - carry out bit(s)

void bitShiftLeft256ymm (__m256i *data, int count)
{
        __m256i innerCarry, rotate;

        innerCarry = _mm256_srli_epi64 (*data, 64 - count);                        // carry outs in bit 0 of each qword
        rotate     = _mm256_permute4x64_epi64 (innerCarry, 0x93);                  // rotate ymm left 64 bits
        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0xFC);   // clear lower qword
        *data    = _mm256_slli_epi64 (*data, count);                               // shift all qwords left
        *data    = _mm256_or_si256 (*data, innerCarry);                            // propagate carrys from low qwords
        //carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        // clear all except lower qword
        //return carryOut;
}

__m256i bitShiftRight256ymm (__m256i *data, int count)
{
        __m256i innerCarry, carryOut, rotate;


        innerCarry = _mm256_slli_epi64(*data, 64 - count);
        rotate =  _mm256_permute4x64_epi64 (innerCarry, 0x39);
        innerCarry = _mm256_blend_epi32 (_mm256_setzero_si256 (), rotate, 0x3F);
        *data = _mm256_srli_epi64(*data, count);
        *data = _mm256_or_si256(*data,  innerCarry);

        carryOut   = _mm256_xor_si256 (innerCarry, rotate);                        //FIXME: not sure if this is correct!!!
        return carryOut;
}
#endif
