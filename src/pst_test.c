
#include <string.h>
#include <zlib.h>

#include "tldevel.h"
#include "tllogsum.h"
#include "tlmisc.h"
#include "tlrng.h"
#include "tlseqio.h"

#include "sim_seq_lib.h"

#define PST_IMPORT
#include "pst.h"

static int test_match_insert(char* infile);

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
        struct kmer_counts* k;
        int i,j,c;
        float out;
        LOG_MSG("%d",argc);
        if(argc == 2){
                filename = argv[1];
                RUN(test_match_insert(filename));
        }else{
                RUN(bar_test());
        }
        //GGTTTACT
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}

int test_match_insert(char* infile)
{
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;
        struct kmer_counts* k = NULL;

        double s[5];
        char* test_seq = NULL;

        float P_M;
        float P_R;
        int test_len = 0;
        int i_point;
        int i,j;
        LOG_MSG("%s",infile);
        if(!my_file_exists(infile)){
                ERROR_MSG("File %s not found");
        }
        RUN(read_10x_white_list(&sb, infile));
        RUN(alloc_kmer_counts(&k, 16));
        RUN(add_counts(k, sb));




        RUN(run_build_pst(&p, k));
        RUN(rm_counts(k,sb));
        RUN(test_kmer_counts(k));
        rng = init_rng(0);

        DECLARE_TIMER(timer);

        for(i = 0; i < 1;i++){

                test_len = 15;

                START_TIMER(timer);
                for(j = 0;j <  5;j++){
                        s[j] = 0.0;
                }
                for (j = 0; j < 1000000;j++){
                        RUN(generate_random_seq(&test_seq, &test_len, rng));
                        RUN(insert_seq(test_seq, test_len, sb->sequences[0]->seq, sb->sequences[0]->len, rng,&i_point));
                        score_pst(p, test_seq, test_len, &P_M,&P_R);

                        s[0]++;
                        s[1] += P_M;
                        s[2] += P_M* P_M;
                }
                s[3] = s[1] / s[0];

                s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
                STOP_TIMER(timer);
                LOG_MSG("Insert: %f %f (in %f sec.)", s[3],s[4], GET_TIMING(timer));

                for(j = 0;j <  5;j++){
                        s[j] = 0.0;
                }
                //START_TIMER(timer);
                for (j = 0; j < 1000000;j++){
                        RUN(generate_random_seq(&test_seq, &test_len, rng));
                        //RUN(insert_seq(test_seq, test_len, sb->sequences[0]->seq, sb->sequences[0]->len, rng));
                        score_pst(p, test_seq, test_len, &P_M,&P_R);

                        s[0]++;
                        s[1] += P_M;
                        s[2] += P_M* P_M;

                }
                s[3] = s[1] / s[0];

                s[4] = sqrt ( (s[0] * s[2] -  pow(s[1], 2.0)) /  (s[0] * ( s[0] - 1.0)));
                START_TIMER(timer);
                for (j = 0; j < 1000000;j++){
                        score_pst(p, test_seq, test_len, &P_M,&P_R);
                }
                STOP_TIMER(timer);
//LOG_MSG("RANDOM: %f %f", s[3],s[4]);
                LOG_MSG("RANDOM: %f %f (in %f sec.)", s[3],s[4], GET_TIMING(timer));


        }

        if(test_seq){
                MFREE(test_seq);
        }
        free_kmer_counts(k);
        //if()
        free_pst(p);
        free_tl_seq_buffer(sb);
        free_rng(rng);
        return OK;
ERROR:
        return FAIL;
}

int bar_test(void)
{
        struct kmer_counts* k = NULL;
        struct tl_seq_buffer* sb = NULL;
        struct pst* p = NULL;
        struct rng_state* rng = NULL;
        double error_rate;
        float P_M,P_R;
        float Q,pbest;
        int i,j,c;
        int addition;

        float back[4];
        float sum;
        char* ref_seq[20] = {
                "AAACCCAAGAAACACT",
                "AAACCCAAGAAACCAT",
                "AAACCCAAGAAACCCA",
                "AAACCCAAGAAACCCG",
                "AAACCCAAGAAACCTG",
                "AAACCCAAGAAACGAA",
                "AAACCCAAGAAACGTC",
                "AAACCCAAGAAACTAC",
                "AAACCCAAGAAACTCA",
                "AAACCCAAGAAACTGC",
                "AAACCCAAGAAACTGT",
                "AAACCCAAGAAAGAAC",
                "AAACCCAAGAAAGACA",
                "AAACCCAAGAAAGCCT",
                "AAACCCAAGAAAGCGA",
                "AAACCCAAGAAAGGAT",
                "AAACCCAAGAAAGGTA",
                "AAACCCAAGAAAGTCT",
                "AAACCCAAGAAATAGG",
                "AAACCCAAGAAATCCA"
        };



        RUN(alloc_tl_seq_buffer(&sb, 16));
        for(i = 0; i < 4;i++){
                back[i] = 0.0f;
        }
        sum = 0.0f;
        for(i = 0; i < 20;i++){
                for(j = 0; j < 16;j++){
                        switch (ref_seq[i][j]) {
                        case 'A': {
                                back[0]++;
                                break;
                        }
                        case 'C': {
                                back[1]++;
                                break;
                        }
                        case 'G': {
                                back[2]++;
                                break;
                        }
                        case 'T': {
                                back[3]++;
                                break;
                        }
                        default:
                                break;
                        }
                        sum++;


                }
        }


        for(i = 0; i < 4;i++){
                back[i] = prob2scaledprob(back[i] / sum);
        }

        //scan_read_with_pst(p,sb->sequences[0]->seq,sb->sequences[0]->len,&out);

        //exit(0);

        RUN(alloc_kmer_counts(&k, 12));
        rng = init_rng(0);


        for(i = 0; i < 20;i++){
                error_rate = (double) i / 20.0f;
                //RUN(score_bar(p, rng, sb->sequences[0]->seq, sb->sequences[0]->len,error_rate));
        }

        struct pst* p_array[20];// = NULL;
        float scores[21];

        for(i = 0; i < 20;i++){
                sb->num_seq = 0;

                snprintf(sb->sequences[sb->num_seq]->seq,sb->sequences[sb->num_seq]->malloc_len,"%s",ref_seq[i]);
                sb->sequences[sb->num_seq]->len = 16;
                sb->num_seq++;

                RUN(add_counts(k,sb));
                RUN(run_build_pst(&p_array[i], k));
                RUN(rm_counts(k,sb));

        }
        //exit(0);
        DECLARE_TIMER(t);
        for(i = 0; i < 20;i++){
                fprintf(stdout,"%s\t",ref_seq[i]);

                score_pst(p_array[0], ref_seq[i], 16, &P_M,&P_R);
                scores[20] = P_R;
                sum = prob2scaledprob(0.0f);
                sum = logsum(sum, P_R);
                START_TIMER(t);
                //for(c = 0; c < 679488;c++){

                        for(j = 0; j < 20;j++){
                                //scan_read_with_pst(p_array[j],ref_seq[i] , 16,&out);
                                //fprintf(stdout,"%f ",out);
                                score_pst(p_array[j], ref_seq[i], 16, &P_M,&P_R);
                                //fprintf(stdout,"(%f)",out);
                                scores[j] = P_M;
                                sum = logsum(sum, P_M);
                        }

                        for(j = 0; j < 21;j++){

                        pbest = 1.0 - scaledprob2prob(scores[j]-sum);

                        if(!pbest){
                                Q = 60.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;
                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }

                        fprintf(stdout,"%3d ",(int) Q);
                }
                fprintf(stdout,"\n");
                        //}
                STOP_TIMER(t);
                //LOG_MSG("Took: %f", GET_TIMING(t));
        }

        char* test_seq = NULL;
        int new_error;
        MMALLOC(test_seq, sizeof(char) * 17);
        for(i = 0; i < 20;i++){
                //int mutate_seq(char* ref,char* target,int len, float error_rate, struct rng_state* rng)
                RUN(mutate_seq(ref_seq[i], test_seq, 16, 0.02, rng,&new_error));
                score_pst(p_array[0], test_seq, 16,  &P_M,&P_R);
                scores[20] = P_R;
                sum = prob2scaledprob(0.0f);
                sum = logsum(sum, P_R);
                fprintf(stdout,"%s\n",ref_seq[i]);
                fprintf(stdout,"%s\t",test_seq);
                START_TIMER(t);
                //for(c = 0; c < 679488;c++){

                for(j = 0; j < 20;j++){
                        //scan_read_with_pst(p_array[j],ref_seq[i] , 16,&out);
                        //fprintf(stdout,"%f ",out);
                        score_pst(p_array[j], test_seq, 16, &P_M,&P_R);
                        //fprintf(stdout,"(%f)",out);
                        scores[j] = P_M;
                        sum = logsum(sum, P_M);
                }
                //fprintf(stdout,"\n");
                for(j = 0; j < 21;j++){

                        pbest = 1.0 - scaledprob2prob(scores[j]-sum);

                        if(!pbest){
                                Q = 60.0;
                        }else if(pbest == 1.0){
                                Q = 0.0;
                        }else{
                                Q = -10.0 * log10(pbest) ;
                        }
                        fprintf(stdout,"%3d ",(int) Q);
                }
                fprintf(stdout,"\n");
                //}
                STOP_TIMER(t);
                //LOG_MSG("Took: %f", GET_TIMING(t));
        }

        MFREE(test_seq);


        for(i = 0; i < 20;i++){
                free_pst(p_array[i]);
        }

        free_rng(rng);

        free_tl_seq_buffer(sb);
        free_kmer_counts(k);
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

