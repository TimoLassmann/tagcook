#include "calibrate_pst.h"


#include "seq_stats.h"
#include "tlrng.h"
#include "tldevel.h"
#include "tllogsum.h"

int calibrate_pst(struct pst* pst, struct sequence_stats_info* si,struct rng_state* rng)
{
        float background[4];
        char* test_seq = NULL;
        int len;

        MMALLOC(test_seq, sizeof(char) * 1024);
        for(i = 0; i < 4;i++){
                background[i] = scaledprob2prob(si->background[i]);
                LOG_MSG("%d %f",i, background[i]);
        }

        for(i = 0; i < 1000000;i++){
                len = (int) tl_random_gaussian(rng, mean_seq_len,stdev_seq_len);

                if(len > 1024){
                        len = 1024;
                }


        }


        return OK;
ERROR:
        return FAIL;
}


