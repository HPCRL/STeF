#ifndef MTTKRP_HARDWIRED_CPP
#define MTTKRP_HARDWIRED_CPP
#include "../inc/mttkrp_hardwired.h"


#ifdef OMP
mutex_array* mutex = NULL;
#endif
int reduce(csf* t, int r, matrix* mat)
{
	memset (mat->val, 0 , sizeof(TYPE)*(mat->dim1)*(mat->dim2));
	int parallel_threshold = 1000;
	if (mat->dim1 > parallel_threshold)
	{
		#pragma omp parallel for
		for(int row_id = 0; row_id < mat->dim1 ; row_id++)
		{
			TYPE* outval = mat->val + row_id * (mat->dim2);
			
			for(int i=0; i<t->num_th ; i++)
			{
				TYPE* reduceval = (t->private_mats[i])->val + row_id * (mat->dim2);
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					outval[y] += reduceval[y];
					//printf("reducing %lf %lf in th %d\n", outval[y], reduceval[y],i);
					//reduceval[y] = 0;
				}
			}	
		}
	}
	else
	{
		for(int row_id = 0; row_id < mat->dim1 ; row_id++)
		{
			TYPE* outval = mat->val + row_id * (mat->dim2);
			
			for(int i=0; i<t->num_th ; i++)
			{
				TYPE* reduceval = (t->private_mats[i])->val + row_id * (mat->dim2);
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					outval[y] += reduceval[y];
					//printf("reducing %lf %lf in th %d\n", outval[y], reduceval[y],i);
					//reduceval[y] = 0;
				}
			}	
		}
	}


	return 0;
}

inline idx_t find_parent(csf* t, idx_t index, int mode)
{
	idx_t start = 0;
	idx_t end = t->mlen[mode-1];

	while(end>start+1) // search within the range [start,end)
	{
		idx_t pivot = (end+start)/2;
		if (t->ptr[mode-1][pivot] > index)
			end = pivot;
		else if (t->ptr[mode-1][pivot+1] < index )
			start = pivot +1;
		else if (t->ptr[mode-1][pivot+1] == index )
		{
			start = pivot + 1;
			end = start + 1;
		}
		else
		{
			start = pivot;
			end = start + 1;	
		}
	}	

	return start;
}


int init_hardwired_2leevl()
{
	#ifdef OMP
	if (mutex == NULL)
	{
		mutex = mutex_alloc_custom((t->mlen)[0] , 16);
		//mutex = mutex_alloc_custom(1024 , 16); // This is what splatt is using
	}
	#endif
	return 0;
}