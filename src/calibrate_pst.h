#ifndef CALIBRATE_PST_H
#define CALIBRATE_PST_H



#include "pst.h"
#include "tlseqio.h"


extern int calibrate_pst(struct pst* pst, struct sequence_stats_info* si,struct rng_state* rng);

#endif
