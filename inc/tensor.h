#ifndef TENSOR_H
#define TENSOR_H

struct tensor_csf
{
	idx_t** ptr;
	idx_t* ptrs;
	idx_t** ind;
	idx_t* inds;
	TYPE* val;
	int* fiber_count;
	int* mlen;
	int nmode;
};

typedef struct tensor_csf csf;

#endif