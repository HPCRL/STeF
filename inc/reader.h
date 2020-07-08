#ifndef READER_H
#define READER_H
#include "../inc/util.h"

struct tensor
{
	idx_t** ptr;
	idx_t* ptrs;
	idx_t** ind;
	idx_t* inds;
	TYPE* val;
};

int read_tensor(const char* file,  struct tensor* res=NULL);

#endif