
#include "tllogsum.h"
#include "core_hmm_functions.h"

#include "hmm.h"
#include "hmm_model_bag.h"

#include "misc.h"

int random_score(struct model_bag* mb,const uint8_t* a,const int len)
{
        int i;
        int c;
        const uint8_t* seqa = a -1;

        float add = prob2scaledprob(1.0 - (1.0 / (float)mb->average_raw_length));
        mb->r_score  = prob2scaledprob(1.0);

        for(i = 1; i <= len;i++){
                c = seqa[i];
                mb->r_score  = mb->r_score  + mb->model[0]->background_nuc_frequency[c] + add;
        }
        mb->r_score  += prob2scaledprob(1.0 / (float)mb->average_raw_length) + mb->r_prior;
        return OK;
}

int backward(struct model_bag* mb,const uint8_t* a,const int len)
{
        int i,j;
        int f,g,c;

        int model_len = 0;

        struct hmm* hmm = 0;
        struct hmm_column* c_hmm_column = 0;
        struct hmm_column* p_hmm_column = 0;

        //float previous_silent[MAX_HMM_SEQ_LEN];


        float* previous_silent = mb->previous_silent;// [MAX_HMM_SEQ_LEN];


        float* psilent;
        float* csilent;


        const uint8_t* seqa = a -1;

        for(j = 0; j < mb->num_models;j++){
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        model_len = mb->model[j]->hmms[f]->num_columns-1;
                        for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
                                c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
                                for(i = 0; i <= len+1;i++){
                                        c_hmm_column->M_backward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->I_backward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->D_backward[i] = prob2scaledprob(0.0);
                                }
                        }
                }
                for(i = 0; i <= len+1;i++){
                        mb->model[j]->silent_backward[i] = prob2scaledprob(0.0f);
                }
        }

        for(i = 0; i <= len+1;i++){
                previous_silent[i] = prob2scaledprob(0.0f);
        }
        previous_silent[len+1] = prob2scaledprob(1.0f);


        mb->model[mb->num_models-1]->silent_backward[len+1] = prob2scaledprob(1.0) + mb->model[mb->num_models-1]->skip;

        for(j = mb->num_models-2 ; j >= 0;j--){
                mb->model[j]->silent_backward[len+1] = mb->model[j+1]->silent_backward[len+1] + mb->model[j]->skip;
        }

        //start with last segment...
        for(j = mb->num_models-1 ; j >= 0;j--){
                if(j == mb->num_models-1 ){
                        psilent = previous_silent;
                }else{
                        psilent = mb->model[j+1]->silent_backward;
                }

                //m_pos-=  mb->model[j]->hmms[0]->num_columns;
                //LOG_MSG("Seg:%d Model%d: POS: %d",j,0,m_pos);


                csilent= mb->model[j]->silent_backward;
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        hmm = mb->model[j]->hmms[f];
                        model_len = mb->model[j]->hmms[f]->num_columns-1;
                        //previous_silent[len+1] = logsum(previous_silent[len+1],  current_silent[len+1]+ mb->model[j] ->skip );
                        //csilent[len+1] =psilent[len+1] + mb->model[j]->skip;
                        //LOG_MSG("Seg:%d Model%d: POS: %d\n",j,f,m_pos);

                        for(i = len ; i > 0;i-- ){
                                //LOG_MSG("Dyn only on: %d %d",s,e);
                                //for(i = len-1 ; i >= 0;i--){ /// DONT MESS WITH THIS - the code takes the first letter to calculate the silent states - EVERYHTHIN is ok..
                                // in the last column an I state will emit a letter not present in the seuqence BUT this will not make it when joined to the forward... messy but works.
                                //	c = (int) seqa[i+1];
                                c = (int)seqa[i+1];
                                c_hmm_column = hmm->hmm_column[model_len];


                                //c_hmm_column->M_backward[i] = psilent[i+1] + mb->model[j]->M_to_silent[f] ;

                                c_hmm_column->M_backward[i] = psilent[i+1] + c_hmm_column->transition[MSKIP];// mb->model[j]->M_to_silent[f] ;


                                //fprintf(stderr," Mback at modellen:%d %f %f\n",i, c_hmm_column->M_backward[i] ,c_hmm_column->transition[MSKIP]);

                                //c_hmm_column->I_backward[i] =  psilent[i+1]+ mb->model[j]->I_to_silent[f] ;

                                c_hmm_column->I_backward[i] =  psilent[i+1] + c_hmm_column->transition[ISKIP];//  mb->model[j]->I_to_silent[f] ;


                                c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->M_backward[i+1] + c_hmm_column->transition[IM] + c_hmm_column->m_emit[c]);

                                c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i] , c_hmm_column->I_backward[i+1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c]);

                                //##################################
                                csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][model_len] + c_hmm_column->m_emit[(int)seqa[i]]);
                                csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][model_len] + c_hmm_column->i_emit[(int)seqa[i]]);
                                //##################################
                                //fprintf(stderr," M:%f I:%f \n", c_hmm_column->M_backward[i] ,c_hmm_column->I_backward[i] );


                                c_hmm_column->D_backward[i] = prob2scaledprob(0.0f);
                                for(g = model_len-1;g >= 0;g--){
                                        c_hmm_column = hmm->hmm_column[g];
                                        //m_pos--;
                                        //LOG_MSG("Seg:%d Model%d: POS: %d",j,f,m_pos);
                                        p_hmm_column = hmm->hmm_column[g+1];

                                        c_hmm_column->M_backward[i]  = p_hmm_column->M_backward[i+1] + p_hmm_column->m_emit[c] + c_hmm_column->transition[MM];


                                        c_hmm_column->M_backward[i] = logsum (c_hmm_column->M_backward[i],psilent[i+1] + c_hmm_column->transition[MSKIP]);

                                        //insert - emit previous symbol etc. etc.
                                        c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i] , c_hmm_column->I_backward[i+1] +c_hmm_column->i_emit[c ]  + c_hmm_column->transition[MI]);

                                        //delete - neex to go to previous columns

                                        c_hmm_column->M_backward[i] = logsum(c_hmm_column->M_backward[i],p_hmm_column->D_backward[i] + c_hmm_column->transition[MD]);

                                        // insert state..
                                        // from previous insertion....
                                        c_hmm_column->I_backward[i] = c_hmm_column->I_backward[i+1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c];

                                        c_hmm_column->I_backward[i] = logsum(c_hmm_column->I_backward[i], psilent[i+1] + c_hmm_column->transition[ISKIP]);
                                        //from previous match state....

                                        c_hmm_column->I_backward[i] = logsum( c_hmm_column->I_backward[i],p_hmm_column->M_backward[i+1] + c_hmm_column->transition[IM] + p_hmm_column->m_emit[c]);
                                        //GRERRRRRRR

                                        //delete state

                                        //from previous delection
                                        c_hmm_column->D_backward[i] = p_hmm_column->D_backward[i] + c_hmm_column->transition[DD];

                                        //from previous match (i.e. gap close

                                        c_hmm_column->D_backward[i] = logsum(c_hmm_column->D_backward[i], p_hmm_column->M_backward[i] + p_hmm_column->m_emit[(int) seqa[i]] + c_hmm_column->transition[DM]);


                                        //##################################
                                        csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][g] + c_hmm_column->m_emit[(int)seqa[i]]);
                                        csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][g] + c_hmm_column->i_emit[(int)seqa[i]]);
                                        //##################################


                                }

                                //##################################

                                //c_hmm_column = hmm->hmm_column[0];
                                // link j+1 to j... dfor silent;
                                //csilent[i] = logsum(csilent[i], c_hmm_column->M_backward[i] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[(int)seqa[i]]);
                                //csilent[i] = logsum(csilent[i], c_hmm_column->I_backward[i]+ mb->model[j]->silent_to_I[f][0] + c_hmm_column->i_emit[(int)seqa[i]]);

                                //##################################

                                //fprintf(stderr,"Looking for Insertyion to silent in segment1: %d	%f\n",f, mb->model[j]->silent_to_I[f]);

                                //this should come from previous state .....
                                csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);

                                //m_pos+=model_len;
                                //LOG_MSG("Seg:%d Model%d: POS: %d\n",j,f,m_pos);
                                //exit(0);

                        }

                }

        }

        mb->b_score = mb->model[0]->silent_backward[1] + mb->m_prior;
        //fprintf(stderr,"SCore:%f	%f\n", mb->b_score , scaledprob2prob(mb->b_score) );

        //fprintf(stderr," BACKWARD:::::::::::\n");

        /*for(j = 0; j < mb->num_models;j++){
          for(i = 0; i <= len;i++){
          fprintf(stderr,"%d	%d	%f\n",j,i,mb->model[j]->silent[i]  );
          }
          }
          exit(0);*/
        /*
          for(j = 0; j < mb->num_models;j++){

          for(f = 0;f < mb->model[j]->num_hmms;f++){
          for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){

          c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
          for(i = 0; i <= len;i++){
					//c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, scaledprob2prob( c_hmm_column->M_backward[i]) , scaledprob2prob( c_hmm_column->I_backward[i]),scaledprob2prob (c_hmm_column->D_backward[i])   );
          }
          }
          }
          }
        */
        //exit(0);
        return OK;
}


struct model_bag* forward(struct model_bag* mb, const uint8_t* a, int len)
{
        int i,j,c;
        int f,g;

        struct hmm* hmm = 0;
        struct hmm_column* c_hmm_column = 0;
        struct hmm_column* p_hmm_column = 0;

        const uint8_t* seqa = a -1;

        float* psilent;
        float* csilent;


        float* previous_silent = mb->previous_silent;

        //float previous_silent[MAX_HMM_SEQ_LEN];
        //float current_silent[MAX_HMM_SEQ_LEN];

        //init

        //float silent_start = prob2scaledprob(1.0);

        // M state of first set of HMMS.....
        for(j = 0; j < mb->num_models;j++){
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
                                c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
                                //i = 0;
                                for(i = 0; i <= len;i++){
                                        c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
                                        //fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, c_hmm_column->M_foward[i] ,c_hmm_column->I_foward[i],c_hmm_column->D_foward[i]    );
                                }
                        }
                }
                for(i = 0; i <= len+1;i++){
                        mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
                }
        }
        //fprintf(stderr,"\n\n\n");
        mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
        //fprintf(stderr,"Init silent states... \n");
        //fprintf(stderr,"%d	%f\n",0,mb->model[0]->silent_forward[0]   );
        for(j = 1; j < mb->num_models;j++){
                mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
                //fprintf(stderr,"%d	%f	%f	%f\n",j,mb->model[j]->silent_forward[0]  ,mb->model[j-1]->silent_forward[0]  , mb->model[j]->skip  );
        }

        for(i = 0; i <= len;i++){
                previous_silent[i] = prob2scaledprob(0.0f);
        }
        previous_silent[0] = prob2scaledprob(1.0);

        //loop thorugh the segments
        // in each run the contained HMMS and update silent states;
        for(j = 0; j < mb->num_models;j++){
                if(j == 0){
                        psilent = previous_silent;
                }else{
                        psilent =  mb->model[j-1]->silent_forward;
                }
                csilent = mb->model[j]->silent_forward;
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        hmm = mb->model[j]->hmms[f];
                        for(i = 1; i <= len;i++){
                                c = seqa[i];

                                c_hmm_column = hmm->hmm_column[0];
                                // first column  comes from previous state cheekily transferring its pd to M[0[
                                c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];

                                c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;

                                //add transitions to first columns////
                                c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);

                                c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];

                                c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);

                                csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);// mb->model[j]->M_to_silent[f]);
                                csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);

                                for(g = 1;g < hmm->num_columns;g++){
                                        c_hmm_column = hmm->hmm_column[g];
                                        p_hmm_column = hmm->hmm_column[g-1];

                                        //Match state

                                        c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];

                                        //transition from previous match state
                                        c_hmm_column->M_foward[i] = logsum( c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
                                        //transition from previous insert state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
                                        //transition from previous delete state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);

                                        // emission promability in curent M state ;
                                        c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];


                                        // Instertion State ..

                                        c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;

                                        //self loop insertion to insertion
                                        c_hmm_column->I_foward[i] =  logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
                                        // start new insertion
                                        c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                        //instertion emission...
                                        c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];
                                        // deletion state
                                        //from previous match state.
                                        c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
                                        //from previous delete state

                                        c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );

                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);// mb->model[j]->M_to_silent[f]);
                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);

                                }
                                //fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);
                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
                                csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
                        }
                }
        }


        mb->f_score = mb->model[mb->num_models-1]->silent_forward[len] + mb->m_prior;

        //fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));

        /*
          for(j = 0; j < mb->num_models;j++){

          for(f = 0;f < mb->model[j]->num_hmms;f++){
          for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){

          c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
          for(i = 0; i <= len;i++){
					//c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
					//c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
					fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i,  scaledprob2prob ( c_hmm_column->M_foward[i]) ,scaledprob2prob ( c_hmm_column->I_foward[i]),scaledprob2prob (c_hmm_column->D_foward[i] )   );
          }
          }
          }
          }*/
        //exit(0);
        return mb;
}



struct model_bag* forward_extract_posteriors(struct model_bag* mb, const uint8_t* a, char* label, int len)
{

        int i,j,c;
        int f,g;

        int hmm_counter = 0;

        struct hmm* hmm = 0;
        struct hmm_column* c_hmm_column = 0;
        struct hmm_column* p_hmm_column = 0;

        const uint8_t* seqa = a -1;

        float* psilent;
        float* csilent;
        float* bsilent;

        //float previous_silent[MAX_HMM_SEQ_LEN];
        //float next_silent[MAX_HMM_SEQ_LEN];
        float* previous_silent = mb->previous_silent;//[MAX_HMM_SEQ_LEN];
        float* next_silent = mb->previous_silent;// [MAX_HMM_SEQ_LEN];


        // M state of first set of HMMS.....
        for(j = 0; j < mb->num_models;j++){
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
                                c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
                                //i = 0;
                                for(i = 0; i <= len;i++){
                                        c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
                                        //fprintf(stderr,"segment:%d	HMM:%d	column:%d	i:%d	%f	%f	%f\n",j,f,g,i, c_hmm_column->M_foward[i] ,c_hmm_column->I_foward[i],c_hmm_column->D_foward[i]    );
                                }
                        }
                }
                for(i = 0; i <= len+1;i++){
                        mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
                }
        }

        //fprintf(stderr,"\n\n\n");
        mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;
        //fprintf(stderr,"Init silent states... \n");
        //fprintf(stderr,"%d	%f\n",0,mb->model[0]->silent_forward[0]   );
        for(j = 1; j < mb->num_models;j++){
                mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
                //mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j-1]->skip ;
        }


        for(i = 0; i <= len;i++){


                next_silent[i] = prob2scaledprob(0.0f);
                previous_silent[i] = prob2scaledprob(0.0f);
        }
        previous_silent[0] = prob2scaledprob(1.0);
        next_silent[len+1] = prob2scaledprob(1.0f);

        //loop thorugh the segments
        // in each run the contained HMMS and update silent states;
        for(j = 0; j < mb->num_models;j++){
                if(j == 0){
                        psilent = previous_silent;
                }else{
                        psilent =  mb->model[j-1]->silent_forward;
                }
                csilent = mb->model[j]->silent_forward;
                if(j +1 != mb->num_models){
                        bsilent = mb->model[j+1]->silent_backward;
                }else{
                        bsilent = next_silent;
                }
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        hmm = mb->model[j]->hmms[f];
                        for(i = 1; i <= len;i++){
                                c = seqa[i];

                                c_hmm_column = hmm->hmm_column[0];
                                // first column  comes from previous state cheekily transferring its pd to M[0[
                                c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];

                                //***************post

                                //if(label[i] == hmm_counter){

                                mb->model[j]->silent_to_M_e[f][0] = logsum(mb->model[j]->silent_to_M_e[f][0] ,psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);


                                c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
                                //}
                                //***************post


                                c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;




                                //add transitions to first columns////
                                c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);

                                c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];

                                //***************post
                                //if(label[i] == hmm_counter){
                                mb->model[j]->silent_to_I_e[f][0]  = logsum(mb->model[j]->silent_to_I_e[f][0] , psilent[i-1] + mb->model[j]->silent_to_I[f][0]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);

                                c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II],c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);

                                c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] +  c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i]-mb->b_score);


                                c_hmm_column->i_emit_e[c] = logsum( c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score);
                                //}
                                //***************post


                                c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);

                                // no post???

                                //

                                csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
                                csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);

                                //***************post
                                //if(label[i] == hmm_counter){
                                c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);

                                c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
                                //}


                                //***************post


                                for(g = 1;g < hmm->num_columns;g++){
                                        c_hmm_column = hmm->hmm_column[g];
                                        p_hmm_column = hmm->hmm_column[g-1];

                                        //Match state


                                        c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];

                                        //transition from previous match state
                                        c_hmm_column->M_foward[i] =  logsum(c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
                                        //transition from previous insert state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
                                        //transition from previous delete state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);

                                        // emission promability in curent M state ;
                                        c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];


                                        //***************post
                                        //if(label[i] == hmm_counter){

                                        //I HOPE THIS IS CORRECT....
                                        mb->model[j]->silent_to_M_e[f][g] = logsum(mb->model[j]->silent_to_M_e[f][g] ,psilent[i-1] + mb->model[j]->silent_to_M[f][g] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);

                                        //I HOPE THIS IS CORRECT....


                                        p_hmm_column->transition_e[MM] = logsum(p_hmm_column->transition_e[MM] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score );

                                        p_hmm_column->transition_e[IM] = logsum(p_hmm_column->transition_e[IM],p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] +  c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);

                                        p_hmm_column->transition_e[DM] = logsum(p_hmm_column->transition_e[DM],p_hmm_column->D_foward[i] + p_hmm_column->transition[DM] + c_hmm_column->m_emit[c] +  c_hmm_column->M_backward[i] -mb->b_score);

                                        c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c],  c_hmm_column->M_foward[i] + c_hmm_column->M_backward[i] -mb->b_score );
                                        //}
                                        //***************post





                                        // Instertion State ..

                                        c_hmm_column->I_foward[i]    =psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;


                                        //self loop insertion to insertion
                                        c_hmm_column->I_foward[i] = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
                                        // start new insertion
                                        c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                        //instertion emission...
                                        c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];

                                        //***************post
                                        //if(label[i] == hmm_counter){
                                        //I HOPE THIS IS CORRECT....

                                        mb->model[j]->silent_to_I_e[f][g]  = logsum(mb->model[j]->silent_to_I_e[f][g] , psilent[i-1] + mb->model[j]->silent_to_I[f][g]  + c_hmm_column->i_emit[c] +  c_hmm_column->I_backward[i] -mb->b_score);

                                        //I HOPE THIS IS CORRECT....


                                        c_hmm_column->transition_e[II] = logsum(c_hmm_column->transition_e[II] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);

                                        c_hmm_column->transition_e[MI] = logsum(c_hmm_column->transition_e[MI] , c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI] + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] - mb->b_score);

                                        c_hmm_column->i_emit_e[c] = logsum(c_hmm_column->i_emit_e[c] , c_hmm_column->I_foward[i]   + c_hmm_column->I_backward[i] - mb->b_score);
                                        //}
                                        //***************post



                                        // deletion state
                                        //from previous match state.
                                        c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
                                        //from previous delete state

                                        c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );

                                        //***************post
                                        //if(label[i] == hmm_counter){
                                        p_hmm_column->transition_e[MD] = logsum(p_hmm_column->transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);

                                        p_hmm_column->transition_e[DD] = logsum(p_hmm_column->transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
                                        //}
                                        //***************post
                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);

                                        //***************post
                                        //if(label[i] == hmm_counter){
                                        c_hmm_column->transition_e[MSKIP] = logsum(c_hmm_column->transition_e[MSKIP], c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP] +  bsilent[i+1] -mb->b_score);

                                        c_hmm_column->transition_e[ISKIP] = logsum(c_hmm_column->transition_e[ISKIP], c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP] +  bsilent[i+1] -mb->b_score);
                                        //}


                                        //***************post



                                }
                                //fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);



                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
                                csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);

                                //***************post
                                //mb->model[j]->M_to_silent_e[f] = logsum(mb->model[j]->M_to_silent_e[f],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f] + bsilent[i+1] -mb->b_score);

                                //fprintf(stderr,"ADDED TO M->S model: %d		%f %f %f %f %f\n", j , scaledprob2prob( mb->model[j]->M_to_silent_e[f]) ,scaledprob2prob(c_hmm_column->M_foward[i]) , scaledprob2prob(mb->model[j]->M_to_silent[f] ),scaledprob2prob( bsilent[i+1]) ,scaledprob2prob( mb->b_score));

                                //mb->model[j]->I_to_silent_e[f] = logsum(mb->model[j]->I_to_silent_e[f] , c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f] + bsilent[i+1] -mb->b_score);
                                //if(label[i] == hmm_counter){
                                mb->model[j]->skip_e =logsum(mb->model[j]->skip_e , psilent[i-1] + mb->model[j]->skip + bsilent[i] -mb->b_score);

                                //}

                                //***************post

                        }
                        hmm_counter++;
                }
        }


        mb->f_score = mb->model[mb->num_models-1]->silent_forward[len] + mb->m_prior;

        //fprintf(stderr,"SCORE:%f	%f\n", mb->f_score, scaledprob2prob(mb->f_score));
        return mb;
}




int forward_max_posterior_decoding(struct model_bag* mb, const uint8_t* a,char** label, int len)
{

        //char* a = ri->seq;
        //int len = ri->len;
        int i,j,c;
        int f,g;

        int hmm_counter = 0;

        struct hmm* hmm = 0;
        struct hmm_column* c_hmm_column = 0;
        struct hmm_column* p_hmm_column = 0;

        const uint8_t* seqa = a -1;

        float* psilent;
        float* csilent;
        //float* bsilent;

        //float previous_silent[MAX_HMM_SEQ_LEN];
        //float next_silent[MAX_HMM_SEQ_LEN];

        float* previous_silent = mb->previous_silent;//[MAX_HMM_SEQ_LEN];
        float* next_silent = mb->previous_silent;// [MAX_HMM_SEQ_LEN];


        ASSERT(mb->b_score != prob2scaledprob(0.0f),"back score is -inf!");
        // M state of first set of HMMS.....
        for(j = 0; j < mb->num_models;j++){
                for(f = 0;f < mb->model[j]->num_hmms;f++){
                        for(g = 0;g < mb->model[j]->hmms[f]->num_columns;g++){
                                c_hmm_column = mb->model[j]->hmms[f]->hmm_column[g];
                                //i = 0;
                                for(i = 0; i <= len;i++){
                                        c_hmm_column->M_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->I_foward[i] = prob2scaledprob(0.0);
                                        c_hmm_column->D_foward[i] = prob2scaledprob(0.0);
                                }
                        }
                }
                for(i = 0; i <= len+1;i++){
                        mb->model[j]->silent_forward[i] = prob2scaledprob(0.0f);
                }
        }

        mb->model[0]->silent_forward[0] = prob2scaledprob(1.0) + mb->model[0]->skip;

        for(j = 1; j < mb->num_models;j++){
                mb->model[j]->silent_forward[0] = mb->model[j-1]->silent_forward[0]  + mb->model[j]->skip ;
        }

        for(i = 0; i <= len;i++){
                for(j = 0; j < mb->total_hmm_num;j++){
                        mb->dyn_prog_matrix[i][j] = prob2scaledprob(0.0f);
                        mb->path[i][j] = -1;
                }
        }


        float total_prob[100];
        for(j = 0; j < mb->total_hmm_num;j++){
                total_prob[j] = prob2scaledprob(0.0);
        }

        for(i = 0; i <= len;i++){
                next_silent[i] = prob2scaledprob(0.0f);
                previous_silent[i] = prob2scaledprob(0.0f);
        }
        previous_silent[0] = prob2scaledprob(1.0);
        next_silent[len+1] = prob2scaledprob(1.0f);


        //loop thorugh the segments
        // in each run the contained HMMS and update silent states;
        for(j = 0; j < mb->num_models;j++){
                if(j == 0){
                        psilent = previous_silent;
                }else{
                        psilent =  mb->model[j-1]->silent_forward;
                }
                csilent = mb->model[j]->silent_forward;
                //	if(j +1 != mb->num_models){
                //		bsilent = mb->model[j+1]->silent_backward;
                //	}else{
                //		bsilent = next_silent;
                //	}
                for(f = 0;f < mb->model[j]->num_hmms;f++){

                        //fprintf(stderr," %d %d %d\n", j , f, hmm_counter);
                        hmm = mb->model[j]->hmms[f];
                        for(i = 1; i <= len;i++){
                                c = seqa[i];

                                c_hmm_column = hmm->hmm_column[0];
                                // first column  comes from previous state cheekily transferring its pd to M[0[
                                c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][0] + c_hmm_column->m_emit[c];

                                //***************post
                                //mb->model[j]->silent_to_M_e[f] = logsum(mb->model[j]->silent_to_M_e[f] ,psilent[i-1] + mb->model[j]->silent_to_M[f] + c_hmm_column->m_emit[c] + c_hmm_column->M_backward[i]  -mb->b_score);
                                //LOG_MSG("%f + %f + %f - %f",total_prob[hmm_counter], c_hmm_column->M_foward[i],c_hmm_column->M_backward[i] ,mb->b_score );
                                total_prob[hmm_counter] = logsum(total_prob[hmm_counter], c_hmm_column->M_foward[i]  +  c_hmm_column->M_backward[i] -mb->b_score );
                                mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );

                                //c_hmm_column->m_emit_e[c] = logsum(c_hmm_column->m_emit_e[c] , c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score);
                                //***************post


                                c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f][0] ;




                                //add transitions to first columns////
                                c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i], c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);

                                c_hmm_column->I_foward[i] = logsum(c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                c_hmm_column->I_foward[i] = c_hmm_column->I_foward[i] + c_hmm_column->i_emit[c];




                                //***************post


                                total_prob[hmm_counter] = logsum(total_prob[hmm_counter], psilent[i-1] + mb->model[j]->silent_to_I[f][0]  + c_hmm_column->i_emit[c] + c_hmm_column->I_backward[i] -mb->b_score );
                                mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );


                                //***************post


                                c_hmm_column->D_foward[i] = prob2scaledprob(0.0f);

                                // no post???

                                //

                                csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + c_hmm_column->transition[MSKIP]);
                                csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + c_hmm_column->transition[ISKIP]);


                                for(g = 1;g < hmm->num_columns;g++){
                                        c_hmm_column = hmm->hmm_column[g];
                                        p_hmm_column = hmm->hmm_column[g-1];

                                        //Match state

                                        //going from silent to internal HMM columns.
                                        c_hmm_column->M_foward[i] = psilent[i-1] + mb->model[j]->silent_to_M[f][g];// + c_hmm_column->m_emit[c];

                                        //transition from previous match state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->M_foward[i-1] + p_hmm_column->transition[MM]);
                                        //transition from previous insert state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i] , p_hmm_column->I_foward[i-1] + p_hmm_column->transition[IM] );
                                        //transition from previous delete state
                                        c_hmm_column->M_foward[i] = logsum(c_hmm_column->M_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DM]);

                                        // emission promability in curent M state ;
                                        c_hmm_column->M_foward[i]  = c_hmm_column->M_foward[i]  + c_hmm_column->m_emit[c];


                                        //***************post
                                        mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->M_foward[i]  + c_hmm_column->M_backward[i] -mb->b_score );
                                        //***************post





                                        // Instertion State ..

                                        c_hmm_column->I_foward[i] = psilent[i-1] + mb->model[j]->silent_to_I[f][g] ;

                                        //self loop insertion to insertion
                                        c_hmm_column->I_foward[i] = logsum (c_hmm_column->I_foward[i] ,  c_hmm_column->I_foward[i-1] + c_hmm_column->transition[II]);
                                        // start new insertion
                                        c_hmm_column->I_foward[i]  = logsum ( c_hmm_column->I_foward[i] ,c_hmm_column->M_foward[i-1] + c_hmm_column->transition[MI]);

                                        //instertion emission...
                                        c_hmm_column->I_foward[i]  = c_hmm_column->I_foward[i]  + c_hmm_column->i_emit[c];

                                        //***************post
                                        mb->dyn_prog_matrix[i][hmm_counter] = logsum(mb->dyn_prog_matrix[i][hmm_counter],  c_hmm_column->I_foward[i] + c_hmm_column->I_backward[i] -mb->b_score );
                                        //***************post



                                        // deletion state
                                        //from previous match state.
                                        c_hmm_column->D_foward[i] = p_hmm_column->M_foward[i] + p_hmm_column->transition[MD];
                                        //from previous delete state

                                        c_hmm_column->D_foward[i] = logsum(c_hmm_column->D_foward[i], p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] );

                                        //***************post
                                        p_hmm_column->transition_e[MD] = logsum(p_hmm_column->transition_e[MD],  p_hmm_column->M_foward[i] + p_hmm_column->transition[MD] + c_hmm_column->D_backward[i] - mb->b_score);

                                        p_hmm_column->transition_e[DD] = logsum(p_hmm_column->transition_e[DD] ,p_hmm_column->D_foward[i] + p_hmm_column->transition[DD] + c_hmm_column->D_backward[i] - mb->b_score);
                                        //***************post

                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] +c_hmm_column->transition[MSKIP]);
                                        csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] +c_hmm_column->transition[ISKIP]);


                                }
                                //fprintf(stderr,"%d	%f	%f	%f	%f	%f\n",i,  c_hmm_column->M_foward[i] , mb->model[j]->M_to_silent[f],  c_hmm_column->I_foward[i] , mb->model[j]->I_to_silent[f], psilent[i]);
                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->M_foward[i] + mb->model[j]->M_to_silent[f]);



                                //csilent[i] =  logsum(csilent[i],  c_hmm_column->I_foward[i] + mb->model[j]->I_to_silent[f]);
                                csilent[i] = logsum(csilent[i], psilent[i] + mb->model[j]->skip);
                        }
                        hmm_counter++;

                }
                //hmm_counter++;
        }

        mb->f_score = mb->model[mb->num_models-1]->silent_forward[len] + mb->m_prior;

        // get barcode score....

        //try to norm ....
        hmm_counter = 0;
        //g = 1;
        next_silent[0] = prob2scaledprob(0.0);
        next_silent[1] = prob2scaledprob(0.0);
        //next_silent[2] = prob2scaledprob(1.0);
        for(j = 0; j < mb->num_models;j++){
                //fprintf(stderr,"MODEL:%d\n",j);
                //next_silent[0] = prob2scaledprob(0.0);
                //next_silent[1] = prob2scaledprob(0.0);

                if(mb->model[j]->num_hmms > 1){
                        g = hmm_counter;
                        next_silent[1] = prob2scaledprob(0.0);


                        for(f = 0;f < mb->model[j]->num_hmms;f++){
                                next_silent[1] = logsum(next_silent[1] , total_prob[hmm_counter]);

                                //fprintf(stderr,"%d %f	%f\n",f,total_prob[hmm_counter],scaledprob2prob( total_prob[hmm_counter]));
                                hmm_counter++;
                        }
                        for(f = 0;f < mb->model[j]->num_hmms;f++){
                                total_prob[g] = total_prob[g]  -next_silent[1];
                                g++;
                        }
                }else{
                        hmm_counter+= mb->model[j]->num_hmms;
                }
        }


        hmm_counter = 0;
        g = 1;
        next_silent[0] = prob2scaledprob(0.0);
        next_silent[1] = prob2scaledprob(0.0);
        next_silent[2] = prob2scaledprob(1.0);

        for(j = 0; j < mb->num_models;j++){
                //fprintf(stderr,"MODEL:%d\n",j);
                //next_silent[0] = prob2scaledprob(0.0);
                //next_silent[1] = prob2scaledprob(0.0);

                if(mb->model[j]->num_hmms > 1){
                        g = 0;
                        next_silent[1] = prob2scaledprob(0.0);


                        for(f = 0;f < mb->model[j]->num_hmms;f++){
                                if(total_prob[hmm_counter] > next_silent[0] && f ){
                                        next_silent[0] = total_prob[hmm_counter];
                                }
                                next_silent[1] = logsum(next_silent[1] , total_prob[hmm_counter]);

                                //fprintf(stderr,"%d %f	%f\n",f,total_prob[hmm_counter],scaledprob2prob( total_prob[hmm_counter]));
                                hmm_counter++;
                        }
                        next_silent[0] = next_silent[0] - next_silent[1]; // this ensures that barprob is never > 1 (happens due to numerical inaccuracy... )
                        //fprintf(stderr,"SUM:%f	%f\n\n", next_silent[1] , scaledprob2prob(next_silent[1] ));

                        next_silent[2] = next_silent[2] +next_silent[0] ;

                }else{
                        hmm_counter+= mb->model[j]->num_hmms;
                }
        }

        if(g){
                mb->bar_score = prob2scaledprob(1.0);
                //ri->bar_prob = prob2scaledprob(1.0);
        }else{
                if(next_silent[2] > 0){
                        mb->bar_score = prob2scaledprob(1.0);
                        //ri->bar_prob = prob2scaledprob(1.0);
                }else{
                        mb->bar_score = next_silent[2];
                        //ri->bar_prob  = next_silent[2];
                }
                //fprintf(stderr,"SELECTED: %f	%f\n", next_silent[2] , scaledprob2prob(next_silent[2] ));
        }

        if(label){
                for(i = 0; i <= len;i++){
                        //	fprintf(stderr,"%d ",i);
                        for(j = 0; j < mb->total_hmm_num;j++){
                                //		fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
                                mb->dyn_prog_matrix[i][j] = scaledprob2prob( mb->dyn_prog_matrix[i][j]);
                                //	total_prob[j] = logsum(total_prob[j] ,  mb->dyn_prog_matrix[i][j]);
                        }

                        //	fprintf(stderr,"\n");
                }
                /*fprintf(stderr,"totalprob: \n");
                  for(j = 0; j < mb->total_hmm_num;j++){
                  fprintf(stderr,"%d	%d	%f	%f\n", j, mb->label[j], total_prob[j], scaledprob2prob(total_prob[j]));
                  }*/


                char* l = *label;
                float max = 0;
                float tmp;
                int move = -1;

                for(i = 1;i <= len;i++){
                        for(j = 0; j < mb->total_hmm_num;j++){
                                max = -1;
                                for(c = 0 ;c <= j ;c++){
                                        tmp =  mb->dyn_prog_matrix[i-1][c] * mb->transition_matrix[c][j];

                                        if(tmp > max){
                                                move = c;
                                                max = tmp;
                                        }
                                        if(tmp == max && c == j){
                                                move = c;
                                                max = tmp;
                                        }

                                        //	fprintf(stderr,"%0.3f ",scaledprob2prob( mb->dyn_prog_matrix[i][j]));
                                }

                                mb->dyn_prog_matrix[i][j]+= max;
                                mb->path[i][j] = move;
                        }
                }
                i = len;
                max = -1;
                for(j = 0; j < mb->total_hmm_num;j++){
                        if(mb->dyn_prog_matrix[i][j] > max){
                                max = mb->dyn_prog_matrix[i][j];
                                move = j;
                        }
                }

                for(i = 0; i <= len;i++){
                        l[i] = 0;
                        //ri->labels[i] = 0;
                }
                l[len] = move;
                //path[len] = move;
                //ri->labels[len] = move;

                for(i = len ;i > 0;i--){
                        move = mb->path[i][move];
                        //	path[i-1] = move;
                        l[i-1] = move;
                        //ri->labels[i-1] = move;
                }
        }
        return OK;
ERROR:
        return FAIL;
}




struct model* reestimate(struct model* m, int mode)
{
        int i,j,c;

        struct hmm_column* m_col = 0;
        //struct hmm_column* copy_col = 0;

        float sum = 0.0f;

        // silent to M /I ....
        // add pseudocount of 1;

        // mode 0
        // train everything...

        //mode 1
        // train everything apart from ssilent to & skip....

        //mode2
        // only train emission probabilities. ....



        if(mode < 1){
                sum = prob2scaledprob(0.0);
                for(i = 0; i < m->num_hmms;i++){
                        for(j = 0; j < m->hmms[i]->num_columns;j++){
                                sum = logsum(sum, logsum(m->silent_to_I_e[i][j],prob2scaledprob(1.0)));
                                sum = logsum(sum, logsum(m->silent_to_M_e[i][j] , prob2scaledprob(1.0)) );
                        }
                        //fprintf(stderr," silent to I: %f",m->silent_to_I_e[i]);
                        //fprintf(stderr," silent to M: %f",m->silent_to_M_e[i]);
                }
                sum = logsum(sum,logsum( m->skip_e , prob2scaledprob(1.0)));
                //fprintf(stderr,"estimated skip: %f\n", m->skip_e);

                for(i = 0; i < m->num_hmms;i++){
                        for(j = 0; j < m->hmms[i]->num_columns;j++){
                                m->silent_to_I[i][j]  =  logsum(m->silent_to_I_e[i][j] ,prob2scaledprob(1.0)) - sum;
                                m->silent_to_M[i][j]  = logsum(m->silent_to_M_e[i][j],prob2scaledprob(1.0)) - sum;
                        }
                }

                m->skip = logsum(m->skip_e ,prob2scaledprob(1.0)) - sum;

                //fprintf(stderr,"SKIP: %f\n", m->skip );

                //clear counts....
                for(i = 0; i < m->num_hmms;i++){
                        for(j = 0; j < m->hmms[i]->num_columns;j++){
                                m->silent_to_I_e[i][j] = prob2scaledprob(0.0);
                                m->silent_to_M_e[i][j] = prob2scaledprob(0.0);
                        }
                }
                m->skip_e = prob2scaledprob(0.0);
        }

        for(i = 0; i < m->num_hmms;i++){
                //copy->silent_to_I[i] = org->silent_to_I[i];
                //copy->silent_to_I_e[i] = org->silent_to_I_e[i];
                //copy->silent_to_M[i] = org->silent_to_M[i];
                //copy->silent_to_M_e[i] = org->silent_to_M_e[i];

                //copy->I_to_silent[i] = org->I_to_silent[i];
                //copy->I_to_silent_e[i] =org->I_to_silent_e[i];
                //copy->M_to_silent[i] = org->M_to_silent[i];
                //copy->M_to_silent_e[i] = org->M_to_silent_e[i];

                for(j = 0; j < m->hmms[0]->num_columns;j++){
                        m_col = m->hmms[i]->hmm_column[j];
                        //copy_col = copy->hmms[i]->hmm_column[j];
                        sum = prob2scaledprob(0.0f);

                        for(c = 0; c< 5;c++){
                                sum = logsum(sum, m_col->i_emit_e[c]);
                        }

                        for(c = 0; c< 5;c++){
                                m_col->i_emit[c] = m_col->i_emit_e[c] - sum;
                        }

                        for(c = 0; c< 5;c++){
                                m_col->i_emit_e[c] =  prob2scaledprob(0.0f);
                        }


                        sum = prob2scaledprob(0.0f);

                        for(c = 0; c< 5;c++){
                                sum = logsum(sum, m_col->m_emit_e[c] );
                        }

                        for(c = 0; c< 5;c++){
                                m_col->m_emit[c] = m_col->m_emit_e[c] - sum;
                        }

                        for(c = 0; c< 5;c++){
                                m_col->m_emit_e[c] =  prob2scaledprob(0.0f);
                        }

                        if(mode < 2){
                                // internal hmm states...
                                if(j != m->hmms[0]->num_columns-1){
                                        sum = prob2scaledprob(0.0f);

                                        //detect &assign
                                        if(m_col->transition[MM] != prob2scaledprob(0.0)){
                                                sum = logsum(sum, m_col->transition_e[MM]);
                                        }
                                        if(m_col->transition[MI] != prob2scaledprob(0.0)){
                                                sum = logsum(sum, m_col->transition_e[MI]);
                                        }
                                        if(m_col->transition[MD] != prob2scaledprob(0.0)){
                                                sum = logsum(sum, m_col->transition_e[MD]);
                                        }
                                        if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
                                                sum = logsum(sum, m_col->transition_e[MSKIP]);
                                        }

                                        //set new prob...
                                        if(m_col->transition[MM] != prob2scaledprob(0.0)){
                                                m_col->transition[MM] = m_col->transition_e[MM]  - sum;
                                        }
                                        if(m_col->transition[MI] != prob2scaledprob(0.0)){
                                                m_col->transition[MI] = m_col->transition_e[MI]  - sum;

                                                //sum = logsum(sum, m_col->transition_e[MI]);
                                        }
                                        if(m_col->transition[MD] != prob2scaledprob(0.0)){
                                                m_col->transition[MD] = m_col->transition_e[MD]  - sum;
                                                //sum = logsum(sum, m_col->transition_e[MD]);
                                        }
                                        if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
                                                m_col->transition[MSKIP] = m_col->transition_e[MSKIP]  - sum;
                                                //sum = logsum(sum, m_col->transition_e[MSKIP]);
                                        }




                                        /*sum = logsum(sum, logsum(m_col->transition_e[MM] , prob2scaledprob(1.0)));
                                          sum = logsum(sum, logsum(m_col->transition_e[MI] , prob2scaledprob(1.0)));
                                          sum = logsum(sum, logsum(m_col->transition_e[MD] , prob2scaledprob(1.0)));
                                          if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
                                          sum = logsum(sum, logsum(m_col->transition_e[MSKIP] , prob2scaledprob(1.0)));
                                          }


                                          m_col->transition[MM] =  logsum(m_col->transition_e[MM] , prob2scaledprob(1.0)) - sum;
                                          m_col->transition[MI] =  logsum(m_col->transition_e[MI] , prob2scaledprob(1.0)) - sum;
                                          m_col->transition[MD] =  logsum(m_col->transition_e[MD], prob2scaledprob(1.0)) - sum;
                                          if(m_col->transition[MSKIP] != prob2scaledprob(0.0)){
                                          m_col->transition[MSKIP] =  logsum(m_col->transition_e[MSKIP], prob2scaledprob(1.0)) - sum;
                                          }
                                        */

                                        sum = prob2scaledprob(0.0f);

                                        sum = logsum(sum, logsum(m_col->transition_e[II] ,prob2scaledprob(1.0)));
                                        sum = logsum(sum, logsum(m_col->transition_e[IM] , prob2scaledprob(1.0)));
                                        if(m_col->transition[ISKIP] != prob2scaledprob(0.0)){
                                                sum = logsum(sum, logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)));
                                        }

                                        m_col->transition[II] =  logsum(m_col->transition_e[II] , prob2scaledprob(1.0)) - sum;
                                        m_col->transition[IM] =  logsum(m_col->transition_e[IM] , prob2scaledprob(1.0)) - sum;
                                        if(m_col->transition[ISKIP] != prob2scaledprob(0.0)){
                                                m_col->transition[ISKIP] =  logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)) - sum;
                                        }
                                        sum = prob2scaledprob(0.0f);

                                        sum = logsum(sum, logsum(m_col->transition_e[DD] , prob2scaledprob(1.0)));
                                        sum = logsum(sum, logsum(m_col->transition_e[DM] , prob2scaledprob(1.0)));

                                        m_col->transition[DD] =  logsum(m_col->transition_e[DD] , prob2scaledprob(1.0)) - sum;
                                        m_col->transition[DM] =  logsum(m_col->transition_e[DM] , prob2scaledprob(1.0)) - sum;



                                }else{ // last hmm column...
                                        // no transitions from M possible....
                                        m_col->transition[MM] =  prob2scaledprob(0.0);
                                        m_col->transition[MI] =  prob2scaledprob(0.0);
                                        m_col->transition[MD] =  prob2scaledprob(0.0);
                                        m_col->transition[MSKIP] = prob2scaledprob(1.0);

                                        //either continue i or goto silent state....
                                        sum = prob2scaledprob(0.0f);

                                        sum = logsum(sum, logsum(m_col->transition_e[II] , prob2scaledprob(1.0)));
                                        sum = logsum(sum, logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)));

                                        //sum = logsum(sum, logsum(m->I_to_silent_e[i] , prob2scaledprob(1.0)));

                                        m_col->transition[II] =  logsum(m_col->transition_e[II] , prob2scaledprob(1.0)) - sum;
                                        m_col->transition[ISKIP] =  logsum(m_col->transition_e[ISKIP] , prob2scaledprob(1.0)) - sum;
                                        //m->I_to_silent[i] =  logsum(m->I_to_silent_e[i] , prob2scaledprob(1.0)) - sum;


                                        //no transtition from D possible.

                                        m_col->transition[DD] = prob2scaledprob(0.0);// m_col->transition_e[DD] + prob2scaledprob(1.0) - sum;
                                        m_col->transition[DM] = prob2scaledprob(0.0);//  m_col->transition_e[DM] + prob2scaledprob(1.0) - sum;


                                }

                                //m->I_to_silent_e[i] = prob2scaledprob(0.0);
                                //m->M_to_silent_e[i] = prob2scaledprob(0.0);

                        }

                        for(c = 0; c< 9;c++){
                                m_col->transition_e[c] =  prob2scaledprob(0.0f);
                        }
                }
        }
        //copy->skip = org->skip;
        //copy->skip_e = copy->skip_e;
        return m;
}



struct hmm* set_hmm_transition_parameters(struct hmm* hmm, int len,double base_error, double indel_freq,  double mean, double stdev)
{
        //cases:
        //1) no MSKIP (apart from last columns
        //2) MSKIP determined by mean & stdev.
        //3 constant MSKIP...


        int i;
        struct hmm_column* col = 0;


        double sum_prob = 0.0;
        if(mean > 0.0 && stdev > 0.0){
                for(i = 0; i <=  len;i++){
                        sum_prob +=gaussian_pdf(i , mean ,stdev);
                }
        }


        if(len == 1){
                //single state - only silent to / from M everything else set to zero....
                col = hmm->hmm_column[0];
                col->transition[MM] = prob2scaledprob(0.0f);
                col->transition[MI] = prob2scaledprob(0.0f);
                col->transition[MD] = prob2scaledprob(0.0f);
                col->transition[MSKIP] = prob2scaledprob(1.0f);

                col->transition[II] = prob2scaledprob(0.0f);
                col->transition[IM] = prob2scaledprob(0.0f);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);
        }else if(len == 2){

                //first column
                col = hmm->hmm_column[0];

                if(mean == -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(0.0);
                }else if(mean > -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
                }else{
                        col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(0 , mean ,stdev) / sum_prob);
                }

                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq ) + prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
                col->transition[MI] = prob2scaledprob(base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));

                col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));


                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);

                //second column
                col = hmm->hmm_column[1];
                col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(1.0);

                col->transition[II] = prob2scaledprob(0.00);
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
                col->transition[DM] = prob2scaledprob(0.0f );//0.999);

        }else{
                //first column....
                col = hmm->hmm_column[0];

                if(mean == -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(0.0);
                }else if(mean > -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
                }else{
                        col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(0 , mean ,stdev) / sum_prob);
                }

                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
                col->transition[MI] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));

                col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));



                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(0.0);

                //middle columns...
                for(i = 1; i < len-2;i++){
                        col = hmm->hmm_column[i];

                        if(mean == -1.0 && stdev == -1.0){
                                col->transition[MSKIP] = prob2scaledprob(0.0);
                        }else if(mean > -1.0 && stdev == -1.0){
                                col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
                        }else{
                                col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf(i , mean ,stdev) / sum_prob);
                        }

                        col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
                        col->transition[MI] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));

                        col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.5)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));



                        col->transition[II] = prob2scaledprob(1.0 - 0.999);
                        col->transition[IM] = prob2scaledprob(0.999);
                        col->transition[ISKIP] = prob2scaledprob(0.0f);


                        col->transition[DD] = prob2scaledprob(1.0 - 0.999);
                        col->transition[DM] = prob2scaledprob(0.999);
                }

                //second last...
                col = hmm->hmm_column[len -2];



                if(mean == -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(0.0);
                }else if(mean > -1.0 && stdev == -1.0){
                        col->transition[MSKIP] = prob2scaledprob(mean / (float)(len-1));
                }else{
                        col->transition[MSKIP] = prob2scaledprob(   gaussian_pdf( len-1.0 , mean ,stdev) / sum_prob);
                }

                col->transition[MM] = prob2scaledprob( 1.0 - base_error * indel_freq)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));
                col->transition[MI] = prob2scaledprob(base_error * indel_freq*1.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));

                col->transition[MD] = prob2scaledprob(base_error * indel_freq*0.0)+ prob2scaledprob (1.0 - scaledprob2prob(col->transition[MSKIP]  ));


                col->transition[II] = prob2scaledprob(1.0 - 0.999);
                col->transition[IM] = prob2scaledprob(0.999);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0);
                col->transition[DM] = prob2scaledprob(1.0);
                //col->transition[DD] = prob2scaledprob(0.0);
                //col->transition[DM] = prob2scaledprob(0.0);

                col = hmm->hmm_column[len -1];

                col->transition[MM] = prob2scaledprob(0.0f );// 1.0 - base_error * indel_freq);
                col->transition[MI] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MD] = prob2scaledprob(0.0f );//(base_error * indel_freq) +  prob2scaledprob(0.5);
                col->transition[MSKIP] = prob2scaledprob(1.0);

                col->transition[II] = prob2scaledprob(0.00);
                col->transition[IM] = prob2scaledprob(0.0);
                col->transition[ISKIP] = prob2scaledprob(0.0f);

                col->transition[DD] = prob2scaledprob(0.0f );//(1.0 - 0.999);
                col->transition[DM] = prob2scaledprob(0.0f );//0.999);
        }
        return hmm;
}




int emit_random_sequence(struct model_bag* mb, uint8_t** r_seq, int* r_len, int t_len, struct rng_state* rng)
{
        uint8_t* s = NULL;
        int s_len;
        //int s_malloc_len = 128;

        int i;
        double r;// = (float)rand()/(float)my_rand_max;
        double sum = prob2scaledprob(0.0f);
        //char alpha[] = "ACGTN";
        int nuc;
        //ri->read_type = 0;

        MMALLOC(s, sizeof(uint8_t) * (t_len+1));
        s_len = 0;

        for(s_len = 0; s_len < t_len;s_len++){
                r = tl_random_double(rng);
                sum = prob2scaledprob(0.0f);
                for(nuc = 0;nuc < 5;nuc++){
                        sum = logsum(sum, mb->model[0]->background_nuc_frequency[nuc] );
                        if(r <  scaledprob2prob(sum)){
                                s[s_len] = nuc;
                                break;
                        }
                }
        }

        s[t_len] = 0;
        //s_len++;

        *r_seq = s;
        *r_len = t_len;

        return OK;
ERROR:
        return FAIL;
}

int emit_read_sequence(struct model_bag* mb, uint8_t** r_seq, int* r_len, int average_length, struct rng_state* rng)
{

        uint8_t* s = NULL;
        int s_len = 0;
        int s_malloc_len = 128;
        int i,j,nuc;
        int state = 0; //0 silent ; 1 M , 2 I , 3 D
        int column = 0;
        int hmm = 0;
        int segment= 0;
        //int status;
        //nt region = 0;
        //int start = 1;
        int len;//mb->model[segment]->hmms[0]->num_columns;
        //char alpha[] = "ACGTN";
        //int parachute = 0;



        double r;

        //fprintf(stderr,"RANd%f MAX:%d\n", r,my_rand_max );

        double sum = prob2scaledprob(0.0f);

        double prob = prob2scaledprob(1.0f);

        MMALLOC(s, sizeof(uint8_t) * s_malloc_len);
        s_len = 0;

        while(s_len <= average_length){
                //KSL_DPRINTF2(("%d %d\n", current_length , average_length ));
                state = 0; //0 silent ; 1 M , 2 I , 3 D
                column = 0;
                hmm = 0;
                segment= 0;

                while(1){

                        //transition
                        r = tl_random_double(rng);
                        //r = (float)rand()/(float)my_rand_max;
                        sum = prob2scaledprob(0.0f);
                        switch (state) {
                        case 0:
                                //	fprintf(stderr,"AM in silent... %f\n",r);
                                len = mb->model[segment]->hmms[0]->num_columns;

                                for(i = 0; i < mb->model[segment]->num_hmms;i++){
                                        for(j = 0;  j < len;j++){
                                                sum = logsum(sum,mb->model[segment]->silent_to_M[i][j]);
                                                //fprintf(stderr,"Trying: hmm:%d Mstate:%d	%f\n", i,j, scaledprob2prob(sum) );
                                                if(r <  scaledprob2prob(sum) ){
                                                        prob += mb->model[segment]->silent_to_M[i][j];
                                                        state = 1;
                                                        column = j;
                                                        hmm = i;
                                                        //if(!start){
                                                        //	segment++;
                                                        //}
                                                        i = 0x7FFFFFF;
                                                        j =  0x7FFFFFF;
                                                        //start = 0;
                                                        break;
                                                }
                                                sum = logsum(sum,mb->model[segment]->silent_to_I[i][j]);
                                                //fprintf(stderr,"Trying: hmm:%d Istate:%d	%f\n", i,j, scaledprob2prob(sum) );
                                                if(r <  scaledprob2prob(sum) ){
                                                        prob += mb->model[segment]->silent_to_I[i][j];
                                                        state = 2;
                                                        column = j;
                                                        hmm = i;
                                                        //if(!start){
                                                        //	segment++;
                                                        //}
                                                        i = 0x7FFFFFF;
                                                        j =  0x7FFFFFF;
                                                        //start = 0;
                                                        break;
                                                }


                                        }
                                }
                                break;
                        case 1:
                                // MM
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MM] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MM];
                                        state = 1;
                                        column++;
                                        //hmm = hmm;
                                        break;
                                }

                                // MI
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MI] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MI];
                                        state = 2;
                                        //column;
                                        //hmm = hmm;
                                        break;
                                }

                                // MD
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MD] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MD];
                                        state = 3;
                                        column++;
                                        //hmm = hmm;
                                        break;
                                }

                                // MSKIP;
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MSKIP] );
                                prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[MSKIP];
                                //if(r <  scaledprob2prob(sum)){
                                state = 0;
                                segment++;
                                column = 0;
                                hmm = 0;
                                //column++;
                                //hmm = hmm;
                                //	break;
                                //}



                                break;

                        case 2:
                                //fprintf(stderr,"PARAM: II %f  IM%f SKIP:%f\n", scaledprob2prob(mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[II] ) ,scaledprob2prob(mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[IM]) , scaledprob2prob(mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[ISKIP] ) );

                                //fprintf(stderr,"%f\n",r);

                                // II
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[II] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[II];
                                        state = 2;
                                        //column++;
                                        //hmm = hmm;
                                        break;
                                }

                                // IM
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[IM] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[IM];
                                        state = 1;
                                        column++;
                                        //hmm = hmm;
                                        break;
                                }

                                // ISKIP
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[ISKIP] );
                                prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[ISKIP];

                                //if(r <  scaledprob2prob(sum)){
                                state = 0;
                                segment++;
                                column = 0;
                                hmm = 0;
                                //column++;
                                //hmm = hmm;
                                //	break;
                                //}


                                break;
                        case 3:
                                // DD
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DD] );

                                if(r <  scaledprob2prob(sum)){
                                        prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DD];
                                        state = 3;
                                        column++;
                                        //hmm = hmm;
                                        break;
                                }

                                // DM
                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DM] );

                                //if(r <  scaledprob2prob(sum)){
                                prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->transition[DM];
                                state = 1;
                                column++;
                                //hmm = hmm;
                                //	break;
                                //}


                                break;


                        default:
                                break;
                        }

                        //fprintf(stderr,"%d segment \n",segment );
                        //fprintf(stderr,"%d hmm \n",hmm );
                        //fprintf(stderr,"%d column \n",column );
                        //fprintf(stderr,"%d state \n",state );
                        //
                        //emit...
                        //r = (float)rand()/(float)my_rand_max;
                        r = tl_random_double(rng);
                        sum = prob2scaledprob(0.0f);

                        if(state == 1){
                                for(nuc = 0;nuc < 5;nuc++){
                                        sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc]);
                                        if(r <  scaledprob2prob(sum)){
                                                prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc];
                                                s[s_len] = nuc;
                                                s_len++;
                                                //fprintf(stderr,"Letter: %d	Segment:%d	hmm:%d	column:%d	state:%d\n",nuc, segment,hmm,column,state );
                                                break;
                                        }
                                }
                                if(nuc == 4){
                                        //fprintf(stderr,"R:%f ",r);
                                        sum = prob2scaledprob(0.0f);
                                        for(nuc = 0;nuc < 5;nuc++){

                                                sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->m_emit[nuc]);
                                                //	fprintf(stderr,"%f ",scaledprob2prob(sum));
                                        }
                                        //fprintf(stderr,"\n");
                                }
                        }

                        if(state == 2){
                                for(nuc = 0;nuc < 5;nuc++){
                                        sum = logsum(sum, mb->model[segment]->hmms[hmm]->hmm_column[column]->i_emit[nuc]);
                                        //		fprintf(stderr,"%d %f %f\n",nuc,r,scaledprob2prob(sum) );
                                        if(r <  scaledprob2prob(sum)){
                                                prob += mb->model[segment]->hmms[hmm]->hmm_column[column]->i_emit[nuc];

                                                s[s_len] = nuc;
                                                s_len++;

                                                break;
                                        }
                                }
                        }

                        if(s_len == s_malloc_len){
                                s_malloc_len = s_malloc_len + s_malloc_len /2;
                                MREALLOC(s, sizeof(uint8_t) * s_malloc_len);
                        }
                        //fprintf(stderr,"segement: %d %d %d\n", segment,mb->num_models,state);
                        if(segment == mb->num_models){
                                break;
                        }
                }
                //fprintf(stderr,"%d len %d \n",current_length,average_length);
                if(s_len >  average_length){
                        //fprintf(stdout,"try again");
                        s_len = 0;
                }else{
                        break;
                }

                //if(current_length+2 >= MAX_HMM_SEQ_LEN){
                //	current_length = 0;
                //}

        }
        s[s_len] = 0;
        //s_len++;

        *r_seq = s;
        *r_len = s_len;
        return OK;
ERROR:
        return FAIL;
}





#ifdef CORE_HMM_TEST

#include "arch_lib.h"

#include "hmm_model_bag.h"

#include "seq_stats.h"

int default_ssi(struct sequence_stats_info** seq_info);
int default_ssi(struct sequence_stats_info** seq_info)
{
        struct sequence_stats_info* ssi = NULL;

        MMALLOC(ssi, sizeof (struct sequence_stats_info));
        ssi->average_length = 500;
        ssi->background[0] = prob2scaledprob(0.25);
        ssi->background[1] = prob2scaledprob(0.25);
        ssi->background[2] = prob2scaledprob(0.25);
        ssi->background[3] = prob2scaledprob(0.25);
        ssi->background[4] = prob2scaledprob(0.0);
        ssi->expected_3_len = NULL;
        ssi->expected_5_len = NULL;
        MMALLOC(ssi->expected_3_len,sizeof(double));
        MMALLOC(ssi->expected_5_len,sizeof(double));
        ssi->expected_3_len[0]  = 0;
        ssi->expected_5_len[0]  = 0;
        ssi->total_num_seq = 1;
        ssi->max_seq_len = 500;
        *seq_info = ssi;
        return OK;
ERROR:
        return FAIL;
}

int main(int argc, char *argv[])
{
        struct arch_library* al = NULL;
        struct model_bag* mb = NULL;
        struct sequence_stats_info* ssi = NULL;
        struct alphabet* a = NULL;
        char* in[] = {
                //"O:N",
                //"S:GTA,AAC",
                "MEREAD:E:N+",
                "I:CCTTAA"
                //"S:ACAGTG,ACTTGA,TTAGGC"
        };
        uint8_t* seq = NULL;
        int seq_len = 8;
        int size;
        int i;

        RUN(create_alphabet(&a, NULL, TLALPHABET_DEFAULT_DNA));
        init_logsum();

        RUN(default_ssi(&ssi));

        RUN(alloc_arch_lib(&al));

        size= sizeof(in) / sizeof(char*);

        read_arch_into_lib(al, in, size);
        LOG_MSG("Read in %d architectures.",al->num_arch);

        MMALLOC(seq, sizeof(uint8_t) * (seq_len+1));
        for(i = 0; i < seq_len;i++){
                seq[i] = i % 4;

        }
        seq[seq_len] =0;

        fflush(stdout);

        RUN(init_model_bag(&mb,al->read_structure[0], ssi,a, 0));




        backward(mb, seq, seq_len);
        forward(mb, seq, seq_len);

        fprintf(stdout,"BScore: %f\n", mb->b_score);
        fprintf(stdout,"FScore: %f\n", mb->b_score);

        free_model_bag(mb);
        MFREE(seq);
        free_arch_lib(al);
        MFREE(ssi->expected_5_len);
        MFREE(ssi->expected_3_len);
        MFREE(ssi);
        free_alphabet(a);
        return EXIT_SUCCESS;
ERROR:
        free_arch_lib(al);
        return EXIT_FAILURE;
}
#endif
