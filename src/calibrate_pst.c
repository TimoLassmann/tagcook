#include "calibrate_pst.h"


#include "seq_stats.h"
#include "tlrng.h"
#include "tldevel.h"
#include "tllogsum.h"
#include "pst_structs.h"
static inline int nuc_to_internal(const char c);

int calibrate_pst(struct pst* pst, struct tl_seq_buffer* sb, int expected_len,struct rng_state* rng)
{

        char* test_seq = NULL;
        double sum;
        double r;
        float P_R,P_M;
        double s0,s1,s2;
        float score;
        int len;
        int i,j,c;

        int p_seq,p_index;
        s0 = 0.0;
        s1 = 0.0;
        s2 = 0.0;

        MMALLOC(test_seq, sizeof(char) * 151);
        for(i = 0; i < 1000000;i++){
                /*p_seq = tl_random_int(rng,sb->num_seq);
                p_index = tl_random_int(rng, sb->sequences[p_seq]->len - expected_len);
                for(j = 0; j < expected_len;j++){
                        test_seq[j] = sb->sequences[p_seq]->seq[p_index+j];
                        }*/
                expected_len = tl_random_int(rng, 100)+50;


                for(j = 0; j < expected_len;j++){
                        sum = 0;
                        r = tl_random_double(rng);
                        for(c = 0; c < 4;c++){
                                sum += scaledprob2prob(pst->fpst_root->prob[0][c]);
                                if(r < sum){
                                        test_seq[j] = "ACGT"[c];
                                        break;
                                }
                        }

                }
                test_seq[expected_len] = 0;
                //LOG_MSG("%s", test_seq);
                RUN(score_pst(pst,test_seq, expected_len, &P_M, &P_R));
                score = P_M-P_R;
                score = score / (double) expected_len;
                //LOG_MSG("%f ",score);
                s0 += 1.0;
                s1 += score;
                s2 += score * score;
        }
        pst->mean = s1 / s0;

        pst->var =  sqrt ( (s0 * s2 -  pow(s1, 2.0)) /  (s0 * ( s0 - 1.0)));
        //exit(0);
        //LOG_MSG("%f %f",pst->mean,pst->var);
        //exit(0);
        MFREE(test_seq);
        return OK;
ERROR:
        return FAIL;
}


static inline int nuc_to_internal(const char c)
{
        switch (c) {
        case 'A':
        case 'a':
                return 0;
                break;
        case 'C':
        case 'c':
                return 1;
                break;
        case 'G':
        case 'g':
                return 2;
                break;
        case 'T':
        case 't':
                return 3;
                break;
        case 'N':
        case 'n':
                return 0;
                break;
        default:
                return 0;
                break;
        }
        return -1;
}
