#ifndef READER_H
#define READER_H
#include "../inc/util.h"

struct tensor_csf
{
	idx_t** ptr;
	idx_t* ptrs;
	idx_t** ind;
	idx_t* inds;
	TYPE* val;
};

typedef struct tensor_csf csf;

int read_tensor(const char* file,  csf* res=NULL);

#endif