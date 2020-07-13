#ifndef TENSOR_H
#define TENSOR_H
#include "../inc/util.h"
struct tensor_csf
{
	idx_t** ptr;
	idx_t* ptrs;
	idx_t** ind;
	idx_t* inds;
	TYPE* val;
	idx_t* fiber_count;
	idx_t* mlen;
	int* modeid;
	int nmode;
	TYPE** intval;
};

struct tensor_coo
{
	idx_t* ind;
	TYPE* val;
	idx_t nnz;
	int nmode;
};

typedef struct tensor_csf csf;
typedef struct tensor_coo coo;


int free_csf(csf* t);
int csf_space(csf* t);
int print_csf(csf* t);


#endif