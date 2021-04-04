#ifndef TENSOR_H
#define TENSOR_H
#include "../inc/util.h"
#include "../inc/matrix.h"
#include "../inc/hash.h"
#include <unordered_set>
#include <array>

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
	idx_t * thread_start;
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



csf* malloc_csf();
coo* malloc_coo();
int free_csf(csf* t);
int free_coo(coo* t);
int csf_space(csf* t);
int print_csf(csf* t,const char* file);
int print_fiber(csf* t, int modeid);
int find_inds(idx_t* inds ,csf* t,idx_t it);
int count_fiber(idx_t** pindex, idx_t nnz, int nmode, int shift, idx_t* fiber_count, int* sort_order);
int coo2csf(coo* dt, csf* t, int* sort_order);
int coo2csf(idx_t** pindex, idx_t* index, TYPE* vals, idx_t nnz, int nmode, idx_t* fiber_count, csf* res,int* mlen, int* sort_order);
int count_fiber(coo* dt, int* sort_order, int hmode);
int count_fiber_leaf_root(idx_t** pindex, idx_t nnz, idx_t modelen, int nmode, int* sort_order);
int count_fiber_leaf_root_fast(csf* t);
#endif