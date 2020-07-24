#ifndef TENSOR_H
#define TENSOR_H
#include "../inc/util.h"
#include "../inc/matrix.h"

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
	matrix** private_mats;
	int num_th;
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
int free_coo(coo* t);
int csf_space(csf* t);
int print_csf(csf* t);
int find_inds(idx_t* inds ,csf* t,idx_t it);


#endif