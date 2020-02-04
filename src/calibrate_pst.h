#ifndef CALIBRATE_PST_H
#define CALIBRATE_PST_H




#include "pst.h"

#include "tlseqio.h"
#include "tlrng.h"

extern int calibrate_pst(struct pst* pst, struct tl_seq_buffer* sb, int expected_len,struct rng_state* rng);

#endif
