
#include <string.h>
#include <zlib.h>

#include "tldevel.h"
#include "tllogsum.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"

#define PST_IMPORT
#include "pst.h"


static int read_10x_white_list(struct tl_seq_buffer** b,char* filename);

static int bar_test(void);
static int score_bar(struct pst* p, struct rng_state*rng, char* ref_seq, int len, double error_rate);

int main(int argc, char *argv[])
{
        char alphabet[] = "ACGT";
        LOG_MSG("Hello World");
        char* filename = NULL;
        char* test_seq = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;
        int i,j,c;
        float out;
        LOG_MSG("%d",argc);
        if(argc == 2){
                filename = argv[1];
                LOG_MSG("%s",filename);
                if(!my_file_exists(filename)){
                        ERROR_MSG("File %s not found");
                }
                RUN(read_10x_white_list(&sb, filename));

                RUN(run_build_pst(&p, sb));

                int num = 0;
                //print_pst(p, p->pst_root, &num);
                LOG_MSG("Found %d leaves ",num);

                num = 0;
                //print_pst(p, p->ppt_root, &num);
                LOG_MSG("Found %d leaves ",num);

                //exit(0);
                for(i = 0; i < 10;i++){
                        fprintf(stdout,">%s\n%s\n",sb->sequences[i]->name,sb->sequences[i]->seq);
                }
                rng = init_rng(0);
                MMALLOC(test_seq, sizeof(char) * (sb->sequences[0]->len+1));
                for(i = 0; i < 10;i++){


                        scan_read_with_pst(p, sb->sequences[i]->seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (overfit)",sb->sequences[i]->seq,out);
                        sb->sequences[i]->seq[6] = 'T';
                        scan_read_with_pst(p, sb->sequences[i]->seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (Added T)",sb->sequences[i]->seq,out);


                        for(j = 0; j < sb->sequences[i]->len;j++){
                                test_seq[j] = alphabet[ tl_random_int(rng,4)];
                        }
                        test_seq[sb->sequences[i]->len]= 0;
                        scan_read_with_pst(p, test_seq, sb->sequences[i]->len,&out);
                        LOG_MSG("%f\t%s (random)",sb->sequences[i]->seq,out);
                }
                //free_error_correct_seq(&e);

                MFREE(test_seq);

                //if()
                free_pst(p);
                free_tl_seq_buffer(sb);
                free_rng(rng);
        }
        RUN(bar_test());
        //GGTTTACT
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}


int bar_test(void)
{

        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;
        double error_rate;
        float out;
        int i,j,c;
        int addition;

        char* ref_seq[10] = {
                "AAACCCAAGAAACACT",
                "AAACCCAAGAAACCAT",
                "AAACCCAAGAAACCCA",
                "AAACCCAAGAAACCCG",
                "AAACCCAAGAAACCTG",
                "AAACCCAAGAAACGAA",
                "AAACCCAAGAAACGTC",
                "AAACCCAAGAAACTAC",
                "AAACCCAAGAAACTCA",
                "AAACCCAAGAAACTGC"
        };

        RUN(alloc_tl_seq_buffer(&sb, 16));
        sb->num_seq = 0;

        snprintf(sb->sequences[sb->num_seq]->seq,sb->sequences[sb->num_seq]->malloc_len,"GGTTTACT");
        sb->sequences[sb->num_seq]->len = 8;
        sb->num_seq++;




        RUN(run_build_pst(&p, sb));


        rng = init_rng(0);

        for(i = 0; i < 100;i++){
                error_rate = (double) i / 100.0f;
                //RUN(score_bar(p, rng, sb->sequences[0]->seq, sb->sequences[0]->len,error_rate));
        }

        struct pst* p_array[10];// = NULL;
        float scores[10];
        float sum;
        for(i = 0; i < 10;i++){
                sb->num_seq = 0;
                snprintf(sb->sequences[sb->num_seq]->seq,sb->sequences[sb->num_seq]->malloc_len,"%s",ref_seq[i]);
                sb->sequences[sb->num_seq]->len = 16;
                sb->num_seq++;
                RUN(run_build_pst(&p_array[i], sb));
        }
        DECLARE_TIMER(t);
        for(i = 0; i < 10;i++){
                sum = prob2scaledprob(0.0f);
                START_TIMER(t);
                //for(c = 0; c < 679488;c++){
                        for(j = 0; j < 10;j++){
                                scan_read_with_pst(p_array[j],ref_seq[i] , 16,&out);
                                //fprintf(stdout,"%f ",out);
                                scores[j] = out;
                                sum = logsum(sum, out);
                        }
                        for(j = 0; j < 10;j++){
                                fprintf(stdout,"%f ", scaledprob2prob(scores[j]-sum));
                        }
                        fprintf(stdout,"\n");
                        //}
                STOP_TIMER(t);
                LOG_MSG("Took: %f", GET_TIMING(t));
        }


        for(i = 0; i < 10;i++){
                free(p_array[i]);
        }


        free_pst(p);
        free_tl_seq_buffer(sb);
        return OK;
ERROR:
        return FAIL;
}

int score_bar(struct pst* p, struct rng_state*rng, char* ref_seq, int len, double error_rate)
{
        char* test_seq = NULL;

        char alphabet[] = "ACGT";
        double s_bar[5];
        double s_time[5];
        int i,j,s;
        int num_tests = 10000;
        float out;
        int num_test_timing = 1000000;
        MMALLOC(test_seq, sizeof(char) * (len+1));
        s_bar[0] = 0.0;
        s_bar[1] = 0.0;
        s_bar[2] = 0.0;

        s_time[0] = 0.0;
        s_time[1] = 0.0;
        s_time[2] = 0.0;
        DECLARE_TIMER(t);
        //START_TIMER(t);
        for(j = 0; j < len;j++){
                out = tl_random_double(rng);
                if(out < error_rate){
                        //LOG_MSG("Has error");
                        test_seq[j] = alphabet[ tl_random_int(rng,4)];
                }else{
                        test_seq[j] = ref_seq[j];
                }
                //test_seq[j] = sb->sequences[0]->seq[j];
        }
        test_seq[len]= 0;

        START_TIMER(t);
        for(j = 0; j < num_test_timing;j++){
                scan_read_with_pst(p, test_seq, len,&out);
        }

        STOP_TIMER(t);
        out = GET_TIMING(t);
        s_time[0]++;
        s_time[1] += out;
        s_time[2] += out * out;

        for(i = 0; i < num_tests;i++){

                //c = tl_random_int(struct rng_state *rng, int a)
                for(j = 0; j < len;j++){
                        out = tl_random_double(rng);
                        if(out < error_rate){
                                //LOG_MSG("Has error");
                                test_seq[j] = alphabet[ tl_random_int(rng,4)];
                        }else{
                                test_seq[j] = ref_seq[j];
                        }
                        //test_seq[j] = sb->sequences[0]->seq[j];
                }
                test_seq[len]= 0;
                ///if(c){
                //LOG_MSG("%f %d %s",out,c,test_seq);

                        ///}
                scan_read_with_pst(p, test_seq, len,&out);
                s_bar[0]++;
                s_bar[1] += out;
                s_bar[2] += out * out;
                //LOG_MSG("%f\t%s (random)",sb->sequences[i]->seq,out);

        }
        s_bar[3] = s_bar[1] / s_bar[0];

        s_bar[4] = sqrt ( (s_bar[0] * s_bar[2] -  pow(s_bar[1], 2.0)) /  (s_bar[0] * ( s_bar[0] - 1.0)));
        s_time[3] = s_time[1] / s_time[0];

        s_time[4] = sqrt ( (s_time[0] * s_time[2] -  pow(s_time[1], 2.0)) /  (s_time[0] * ( s_time[0] - 1.0)));

        LOG_MSG("%f %f  at error %f in %f (+/- %f)", s_bar[3],s_bar[4], error_rate, s_time[3],s_time[4]);
        MFREE(test_seq);
        return OK;
ERROR:
        return FAIL;
}

int read_10x_white_list(struct tl_seq_buffer** b,char* filename)
{
        struct tl_seq_buffer* sb = NULL;
        struct tl_seq* s = NULL;
        gzFile f_ptr;
        char* buffer = NULL;
        char* tmp = NULL;
        int buffer_len = 256;

        MMALLOC(buffer, sizeof(char) * buffer_len);

        RUN(alloc_tl_seq_buffer(&sb, 1000000));


        RUNP(f_ptr = gzopen(filename, "r"));
        while((tmp =  gzgets(f_ptr, buffer, buffer_len)) != NULL){
                //fprintf(stdout,"%s",buffer);
                s = sb->sequences[sb->num_seq];
                s->len = strnlen(buffer, buffer_len) -1;

                snprintf(s->name, TL_SEQ_MAX_NAME_LEN, "Bar%d", sb->num_seq+1);

                while(s->malloc_len < s->len){
                        RUN(resize_tl_seq(s));
                }
                strncpy(s->seq, buffer, s->len);
                s->seq[s->len] = 0;
                //snprintf(s->seq, s->malloc_len, "%s",buffer);
                sb->num_seq++;
                //if(sb->num_seq == 100000){
                //break;
                //}
                if(sb->num_seq == sb->malloc_num){
                        RUN(resize_tl_seq_buffer(sb));
                }
        }

        *b = sb;

        gzclose(f_ptr);


        MFREE(buffer);
        return OK;

ERROR:
        if(buffer){
                MFREE(buffer);
        }
        if(f_ptr){
                gzclose(f_ptr);
        }
        return FAIL;

}
