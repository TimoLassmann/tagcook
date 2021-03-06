#ifndef INIT_HMM_H
#define INIT_HMM_H


#include "hmm.h"

#include "tldevel.h"

#include "tlalphabet.h"
#include "arch_lib.h"

/* generic initialisation */
extern int set_emission_p(struct model* m, const struct segment_specs* spec, const struct alphabet* a, const double* background, const double base_error, const double indel_freq);
extern int set_transition_p(struct model* m, const double base_error, const double indel_freq);

/* specific initialization */
extern int init_plus_model(struct model* m, const int expected_len);


extern struct model* init_model_according_to_read_structure(struct model* model,struct read_structure* rs ,struct alphabet*a,  int key,const  double* background,int assumed_length);

extern struct model* malloc_model_according_to_read_structure(int num_hmm, int length,int dyn_length);


extern void free_model(struct model* model);
#endif
