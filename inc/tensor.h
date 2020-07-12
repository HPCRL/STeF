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
	int* fiber_count;
	int* mlen;
	int* modeid;
	int nmode;
	TYPE** intval;
};

typedef struct tensor_csf csf;


int free_csf(csf* t);
int csf_space(csf* t);

#endif