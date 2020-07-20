#ifndef TENSOR_CPP
#define TENSOR_CPP


#include "../inc/tensor.h"

int free_csf(csf* t)
{
	int i;
	rem(t->ptr);
	rem(t->ptrs);
	rem(t->ind);
	rem(t->inds);
	rem(t->val);
	rem(t->fiber_count);
	rem(t->mlen);
	rem(t->modeid);
	//printf("%d\n", t->nmode );
	if(t->intval != NULL)
		for (i=0; i< (t->nmode) ; i++)
			rem(t->intval[i]);
	rem(t->intval);
	rem(t);
	return 0;
}

int free_coo(coo* t)
{
	rem(t->ind);
	rem(t->val);
	return 0;
}


int print_csf(csf* t)
{
	for(int i=0; i<t->nmode ;i++)
	{
		printf("Mode %d fiber count %d\n", i, t->fiber_count[i]);
		printf("Mode %d mode length %d\n", i, t->mlen[i]);
		//printf("Mode %d fiber_count %d\n", fiber_count[i]);
	}
	return 0;
}

int csf_space(csf* t)
{

	return 0;
}

#endif