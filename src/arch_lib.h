#ifndef ARCH_LIB_H
#define ARCH_LIB_H

#include "pst.h"
#include "correct.h"


#include <stdint.h>

#define ARCH_ETYPE_EXTRACT 1
#define ARCH_ETYPE_APPEND 2
#define ARCH_ETYPE_SPLIT 3
#define ARCH_ETYPE_IGNORE 4
#define ARCH_ETYPE_PARTIAL 5
#define ARCH_ETYPE_WOBBLE_LEFT 6
#define ARCH_ETYPE_WOBBLE_RIGHT 7

#define ARCH_ETYPE_APPEND_CORRECT 8
#define ARCH_ETYPE_CORRECT 9



struct segment_specs{
        khash_t(exact)* bar_hash;
        struct pst* pst;
        char* name;
        char* correct_name;
        char* qual_name;
        char** seq;
        int num_seq;
        int max_len;
        int min_len;
        int alloc_len;
        uint8_t extract;
};

struct read_structure{
        struct segment_specs** seg_spec;
        int num_segments;
        int alloc_num_seqments;
};

struct hash_store{
        khash_t(exact)** hash;
        struct pst** pst;
        char** name;
        int* len;
        int alloc_hash;
        int num_hash;
};

struct arch_library{
        struct read_structure** read_structure;
        char* name;
        char** spec_line;
        float** arch_posteriors;
        float* confidence_thresholds;
        int* arch_to_read_assignment;
        float P;
        uint8_t read_order_check;
        int priority;
        int num_arch;
        int alloc_num_arch;
        int num_file;
};

struct cookbook{
        struct arch_library** lib;
        struct hash_store* hs;
        float* scores;
        int num_lib;
        int alloc_num_lib;
        int best;
};


extern int read_cookbook_command_line(struct cookbook** cookbook, char* in);
extern int read_cookbook_file(struct cookbook** cookbook, char* filename);

extern int free_cookbook(struct cookbook** cookbook);

extern int read_architecture_files(struct arch_library* al, char* filename);
extern int read_arch_into_lib(struct arch_library* al, char** list, int len);

extern int alloc_arch_lib(struct arch_library** arch);
extern void free_arch_lib(struct arch_library* arch);

extern int print_segment_spec(const struct segment_specs* spec);
/* emit sequences directly from read structure  */




#endif
