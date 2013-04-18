//
//  barcode_hmm.h
//  tagdust2
//
//  Created by lassmann on 2/5/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//




#ifndef _MM_ALIGN16
#ifdef __GNUC__
#define _MM_ALIGN16 __attribute__((aligned (16)))
#endif
#ifdef __MSVC__
#define _MM_ALIGN16 __declspec(align(16))
#endif
#endif




#ifndef tagdust2_barcode_hmm_h
#define tagdust2_barcode_hmm_h

#define MM 0
#define MI 1
#define MD 2
#define MQUIT 7
#define II 3
#define IM  4
#define DD 5
#define DM 6


#define SELF 0
#define NEXT 1


#define MAX_HMM_SEQ_LEN 150

#define MAX_NUM_SUB_MODELS 64


struct hmm_column{
	float M_foward[MAX_HMM_SEQ_LEN];
	float M_backward[MAX_HMM_SEQ_LEN];
	
	float I_foward[MAX_HMM_SEQ_LEN];
	float I_backward[MAX_HMM_SEQ_LEN];
	
	float D_foward[MAX_HMM_SEQ_LEN];
	float D_backward[MAX_HMM_SEQ_LEN];
	
	float long_transition[MAX_NUM_SUB_MODELS];
	float long_transition_e[MAX_NUM_SUB_MODELS];
	
	float short_transition[8];
	float short_transition_e[8];
	
	float m_emit[5];
	float i_emit[5];
	float m_emit_e[5];
	float i_emit_e[5];
	
	int identifier;
}_MM_ALIGN16;

struct hmm{
	struct hmm_column** hmm_column;
	int num_columns;
	
}_MM_ALIGN16;

struct model{
	struct hmm** hmms;
	float background_nuc_frequency[5];
	float* silent_to_M;
	
	float* M_to_silent;
	
	float* silent_to_I;
	
	float* I_to_silent;
	
	
	float silent[MAX_HMM_SEQ_LEN];
	float random_next;
	float random_self;
	float skip;
	
	int average_length;
	int num_hmms;
}_MM_ALIGN16;

struct model_bag{
	struct model** model;
	int num_models;
}_MM_ALIGN16;



void hmm_controller(struct parameters* param,int (*fp)(struct read_info** ,struct parameters*,FILE* ),int file_num);

struct model* malloc_model(int main_length, int sub_length, int number_sub_models);

struct model* init_model(struct model* model);
struct model* copy_and_malloc_model(struct model* org);
struct model* add_estimates_to_model(struct model* target, struct model* source);
void free_model(struct model* model);

struct model* malloc_model_according_to_read_structure(struct read_structure* rs, int key);
struct model* init_model_according_to_read_structure(struct model* model,struct parameters* param , int key, float* background,int assumed_length);
void print_model(struct model* model);


struct model_bag* forward(struct model_bag* mb, char* a, int len);
struct model_bag* backward (struct model_bag* mb, char* a, int len);

#endif




