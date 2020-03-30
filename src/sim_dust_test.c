#include <getopt.h>
#include "tldevel.h"
#include "tlseqio.h"

#include "tlrng.h"
#include "sim_seq_lib.h"
#include <string.h>

#define OPT_ERROR_RATE 1
#define OPT_INDEL_FRAC 2


static void usage();

struct param{
        char* ref_name;
        char* rRNA_name;
        char* out_name;
        char* mode;
        int seed;
        double error_rate;
        double indel_frac;
        int n_sim_seq;
        int sim_len;

};

int init_param(struct param**p);
void free_param(struct param* p);
int eval(struct param* p);

static int sim_reads(struct param* p);

int main(int argc, char *argv[])
{
        struct param* p = NULL;
        int c;
        if (argc < 2){
                usage();
                return OK;
        }

        RUN(init_param(&p));

        while (1){

                static struct option long_options[] ={
                        {"erate", required_argument, 0, OPT_ERROR_RATE},
                        {"indelfrac", required_argument, 0, OPT_INDEL_FRAC},
                        {0, 0, 0, 0}
                };

                int option_index = 0;
                c = getopt_long_only (argc, argv,"m:t:i:r:o:s:",long_options, &option_index);

                if (c == -1){
                        break;
                }

                switch(c) {
                case 0:
                        break;
                case OPT_ERROR_RATE:
                        p->error_rate = atof (optarg);
                        break;
                case OPT_INDEL_FRAC:
                        p->indel_frac = atof (optarg);
                        break;
                case 'm':
                        p->mode = optarg;
                        break;
                case 't':
                case 'i':
                        p->ref_name = optarg;
                        break;
                case 'o':
                        p->out_name = optarg;
                        break;
                case 'r':
                        p->rRNA_name = optarg;
                        break;
                case 's':
                        p->seed= atoi(optarg);
                        break;
                case '?':
                        exit(1);
                        break;
                default:
                        ERROR_MSG("Option not recognised");
                }
        }

        if(!p->mode){
                ERROR_MSG("No mode specified");

        }

        if(!strncmp(p->mode,"sim",3)){

                if(!p->ref_name){
                        usage();
                        ERROR_MSG("no ref name");
                }

                if(!p->rRNA_name){
                        usage();
                        ERROR_MSG("no rRNA name");
                }

                if(!p->out_name){
                        usage();
                        ERROR_MSG("no out name");
                }
                //MMALLOC(p, size)

                RUN(sim_reads(p));
        }else if(!strncmp(p->mode,"eval", 4)){
                LOG_MSG("Eval mode");
                RUN(eval(p));
        }else{
                ERROR_MSG("Mode: %s not recognised",p->mode);
        }

        free_param(p);
        return EXIT_SUCCESS;
ERROR:
        free_param(p);
        return EXIT_FAILURE;
}

int eval(struct param* p)
{
        struct tl_seq_buffer* r = NULL;
        struct file_handler* r_file = NULL;
        char buffer[128];
        int class, n_t, n_f;
        int i = 0;
        double tp,fp,tn,fn;

        RUN(open_fasta_fastq_file(&r_file,p->ref_name, TLSEQIO_READ));
        RUN(read_fasta_fastq_file(r_file, &r, 1000000));
        LOG_MSG("Read %d sequences from %s",r->num_seq, p->ref_name);


        close_seq_file(&r_file);
        if(!r->num_seq){
                free_tl_seq_buffer(r);
                ERROR_MSG("No sequences found!");
        }
        tp = 0.0;
        fp = 0.0;
        tn = 0.0;
        fn = 0.0;
        n_t = 0;
        n_f = 0;
        class = -1;
        sscanf(r->sequences[0]->name, "%*s %d %d %d",&class, &n_t,&n_f);

        tn = (double) n_t;
        fn = (double) n_f;


        for(i = 0; i < r->num_seq;i++){

                sscanf(r->sequences[i]->name, "%*s %d %d %d",&class, &n_t,&n_f);
                //snprintf(s_buffer->sequences[c]->name,TL_SEQ_MAX_NAME_LEN,"%s_%d_%d_%d","gene",1, p->n_sim_seq/2,p->n_sim_seq/2);
                //fprintf(stdout,"%d %s %d %d %d\n",i,r->sequences[i]->name ,class,n_t,n_f);
                if(class == 1){
                        //fn--;
                        tp++;
                }else if(class == 0){
                        fp++;
                        //tn--;
                }


        }
        fn = n_t - tp;
        tn = n_f - fp;
        free_tl_seq_buffer(r);
        LOG_MSG("TP: %f", tp);
        LOG_MSG("TN: %f", tn);
        LOG_MSG("FP: %f", fp);
        LOG_MSG("FN: %f", fn);
        return OK;
ERROR:
        return FAIL;
}

int sim_reads(struct param* p)
{
        struct tl_seq_buffer* r = NULL;
        struct tl_seq_buffer* b = NULL;
        struct tl_seq_buffer* s_buffer = NULL;
        struct file_handler* r_file = NULL;
        struct file_handler* b_file = NULL;
        struct file_handler* o_file = NULL;
        struct rng_state* rng = NULL;
        char* tmp_seq = NULL;

        int i,j,c;
        int s_index;
        int pos;
        int errors;

        RUNP(rng = init_rng(p->seed));

        RUN(open_fasta_fastq_file(&r_file,p->ref_name, TLSEQIO_READ));
        RUN(read_fasta_fastq_file(r_file, &r, 1000000));
        LOG_MSG("Read %d sequences from %s",r->num_seq, p->ref_name);

        close_seq_file(&r_file);


        RUN(open_fasta_fastq_file(&b_file,p->rRNA_name, TLSEQIO_READ));
        RUN(read_fasta_fastq_file(b_file, &b, 1000000));
        LOG_MSG("Read %d sequences from %s",b->num_seq, p->rRNA_name);
        close_seq_file(&b_file);

        RUN(alloc_tl_seq_buffer(&s_buffer, p->n_sim_seq));
        c = 0;
        for(i = 0; i < p->n_sim_seq/2;i++){
                s_index = tl_random_int(rng, r->num_seq);
                while(r->sequences[s_index]->len < p->sim_len){
                        s_index = tl_random_int(rng, r->num_seq);
                }
                pos =  tl_random_int(rng, r->sequences[s_index]->len - p->sim_len);

                while(p->sim_len+1 > s_buffer->sequences[c]->malloc_len){
                        resize_tl_seq(s_buffer->sequences[c]);
                }
                for(j = 0; j < p->sim_len;j++){
                        s_buffer->sequences[c]->seq[j] = r->sequences[s_index]->seq[j+pos];
                        s_buffer->sequences[c]->qual[j] = 'A';
                }
                s_buffer->sequences[c]->seq[p->sim_len] = 0;
                s_buffer->sequences[c]->qual[p->sim_len] = 0;
                snprintf(s_buffer->sequences[c]->name,TL_SEQ_MAX_NAME_LEN,"%s %d %d %d","gene",1, p->n_sim_seq/2,p->n_sim_seq/2);
                s_buffer->sequences[c]->len = p->sim_len;
                //LOG_MSG("%s",s_buffer->sequences[c]->seq);
                c++;

        }
        for(i = 0; i < p->n_sim_seq/2;i++){
                s_index = tl_random_int(rng, b->num_seq);
                while(b->sequences[s_index]->len < p->sim_len){
                        s_index = tl_random_int(rng, b->num_seq);
                }
                pos =  tl_random_int(rng, b->sequences[s_index]->len - p->sim_len);

                while(p->sim_len+1 > s_buffer->sequences[c]->malloc_len){
                        resize_tl_seq(s_buffer->sequences[c]);
                }
                for(j = 0; j < p->sim_len;j++){
                        s_buffer->sequences[c]->seq[j] = b->sequences[s_index]->seq[j+pos];
                        s_buffer->sequences[c]->qual[j] = 'A';
                }
                s_buffer->sequences[c]->seq[p->sim_len] = 0;
                s_buffer->sequences[c]->qual[p->sim_len] = 0;
                snprintf(s_buffer->sequences[c]->name,TL_SEQ_MAX_NAME_LEN,"%s %d %d %d","rRNA",0, p->n_sim_seq/2,p->n_sim_seq/2);
                s_buffer->sequences[c]->len = p->sim_len;
                //snprintf(s_buffer->sequences[c]->name,TL_SEQ_MAX_NAME_LEN,"%s%d","rRNA",c);
                //LOG_MSG("%s (%d)",s_buffer->sequences[c]->seq, s_buffer->sequences[c]->len);
                c++;
        }
        s_buffer->num_seq = c;
        s_buffer->is_fastq = 1;

        /* mutate */
        MMALLOC( tmp_seq,sizeof(char) * (1+ p->sim_len));
        for(i = 0; i < s_buffer->num_seq;i++){


                RUN(mutate_seq(s_buffer->sequences[i]->seq, tmp_seq, s_buffer->sequences[i]->len, p->error_rate, rng, &errors));
                LOG_MSG("errors: %d (%f)",errors, p->error_rate);
                for(j = 0; j < s_buffer->sequences[i]->len;j++){
                        s_buffer->sequences[i]->seq[j] = tmp_seq[j];
                }
        }

        RUN(open_fasta_fastq_file(&o_file, p->out_name, TLSEQIO_WRITE));
        write_seq_buf(s_buffer, o_file);
        close_seq_file(&o_file);

        free_tl_seq_buffer(s_buffer);
        free_tl_seq_buffer(r);
        free_tl_seq_buffer(b);

        free_rng(rng);
        return OK;
ERROR:
        return FAIL;
}

void usage()
{
        fprintf(stdout, "Usage:   sim_dust -m <sim/eval> -t <gencode.fasta> -r <ribosomal ref> -o outfilename \n\n");
        fprintf(stdout, "Options:\n");


        fprintf (stdout,"\t%-17s%10s%7s%-30s\n","-Q","FLT","", "confidence threshold [20].");
}


int init_param(struct param**p)
{
        struct param* param = NULL;
        MMALLOC(param, sizeof(struct param));
        param->error_rate = 0.02;
        param->indel_frac = 0.1;
        param->seed = 0;
        param->out_name = NULL;
        param->rRNA_name = NULL;
        param->ref_name = NULL;
        param->n_sim_seq = 100000;
        param->sim_len = 75;
        *p = param;
        return OK;
ERROR:
        return FAIL;
}

void free_param(struct param* p)
{
        if(p){
                MFREE(p);
        }
}
