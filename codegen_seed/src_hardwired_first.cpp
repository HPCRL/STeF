#ifndef MTTKRP_HARDWIRED_CPP
#define MTTKRP_HARDWIRED_CPP
#include "../inc/mttkrp_hardwired.h"

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
