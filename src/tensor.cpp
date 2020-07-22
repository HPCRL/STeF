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

int find_inds(idx_t* inds ,csf* t,idx_t it)
{
	int nmode = t->nmode;
	//printf("here ss\n");
	if (it == 0)
	{
		for (int i=0; i<nmode ; i++)
		{
			inds[i] = 0;
		}

		return 0;
	}
	else
	{
		
		inds[nmode-1] = it;
		idx_t last_pos = it;
		for(int i=nmode-2 ; i>=0 ; i--)
		{
			// do binary seach and find the index
			// idx_t id = -1;
			idx_t start = 0;
			idx_t end = t->fiber_count[i];
			while(end > start+1)
			{
				idx_t pivot = (end + start)/2;
				if(t->ptr[i][pivot] > last_pos)
				{
					end = pivot;
				}
				else if (t->ptr[i][pivot] <= last_pos)
				{
					start = pivot;
				}
				//printf("%d %d %d %d\n",start, end , i , last_pos);
			}
			inds[i] = start;
			last_pos = start;
		}
		/*
		for(int i=0; i<nmode ; i++)
		{
			printf("%d", inds[i] );
		}
		printf("\n");
		*/
		return 0;
	}

	
}


#endif