#ifndef MTTKRP_CPP
#define MTTKRP_CPP
#include "../inc/mttkrp.h"



int mttkrp_atomic3(csf* t, int mode, int r, matrix** mats)
{
	TYPE* partial_products;
	int i,j,k,x,y,ind,nmode;
	TYPE* vals;
	partial_products = (TYPE* ) malloc(nmode*r*sizeof(TYPE));
	vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	nmode = t->nmode;
	matrix m;

	for(i=0 ; i<mats[mode]->dim1*mats[mode]->dim2 ; i++)
		vals[i] = 0;

	for(i = 0 ; i < t->fiber_count[0] ; i++ )
	{
		
		for(y = 0; y<r ; y++)
		{
			partial_products[y] = MAT(mats[0], t->ind[0][i] ,y);
		}
		for(j = t->ptr[0][i] ; j < t->ptr[0][i+1] ; j++)
		{

			for(y = 0; y<r ; y++)
			{
				partial_products[r + y] = partial_products[y] * MAT(mats[1], t->ind[1][j] ,y);
			}

			for(k = t->ptr[1][j] ; k < t->ptr[1][j+1]; k++)
			{
				int kk = t->ind[2][k];
				for(y = 0; y<r ; y++)
				{
					vals[kk*mats[mode]->dim2 + y] += partial_products[r + y] * t->val[k];
					if (y == 0)
						printf("%lf %lf %d %d\n", partial_products[r + y] * t->val[k] , vals[kk*mats[mode]->dim2 + y], kk*mats[mode]->dim2 + y, kk);
					//printf("%lf\n", t->val[k]);
				}
			}
		}
		
	}

	mats[nmode-1]->val = vals;
	return 0;
}


int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats)
{
	/* 
	if(t->nmode == 3)
	{
		return mttkrp_atomic3(t,mode,r,mats);
	}
	*/

	

	TYPE* partial_products, *vals;
	int i,ii,it,nmode,nnz;
	idx_t* inds;
	int num_th, th;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
	#endif

	printf("num ths %d\n", num_th);
	partial_products = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	inds = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	nnz = t->fiber_count[nmode-1];
	vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	for(i=0 ; i<mats[mode]->dim1*mats[mode]->dim2 ; i++)
		vals[i] = 0;
	


	it = 0;
	// DO the process for the first nnz

	for(i = 0 ; i<r ; i++)
	{
		partial_products[i] = MAT(mats[0], t->ind[0][0] ,i);
	}

	for(ii = 1; ii < nmode - 1 ; ii++)
	{
		for(i = 0 ; i<r ; i++)
		{
			partial_products[ii * r + i]  = partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][0],i);
		}	
	}

	for(i=0 ; i<r ; i++)
	{
		vals[t->ind[nmode-1][0]*r + i]	= partial_products[(nmode-2) * r + i] * t->val[0];
	}

	for(i = 0 ; i < nmode ; i++)
	{
		inds[i] = 0;
	}

	it = 1;
	printf("%d\n",nnz);
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		#ifdef OMP
			th = omp_get_thread_num();
		#endif
		while(it<nnz)
		{
			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				inds[th*nmode + nmode-2] ++;
				for(i = nmode-3; i>=0 ; i--)
				{
					if(inds[th*nmode + i+1]  < t->ptr[i][inds[th*nmode + i]+1] )
					{
						update = i+1;
						break;
					}
					else
					{
						inds[i] ++;
	
					}
	
				}
				if(update == 0)
				{
					for(i = 0 ; i<r ; i++)
					{
						partial_products[th*nmode*r + i] = MAT(mats[0], t->ind[0][inds[nmode*th  +0]] ,i);
					}
					update ++;
				}
	
	
	
				for(ii = update; ii<nmode-2; ii++)
				{
					for(i = 0 ; i<r ; i++)
					{
						partial_products[th*nmode*r +  ii*r + i] = MAT(mats[ii], t->ind[ii][inds[th*nmode + ii]] ,i);
					}
					update ++;
				}
			}
			// use partial products to update last mode
			for(i=0 ; i<r ; i++)
			{
				// put a locking step here
				// This should be atomic
				vals[t->ind[nmode-1][it]*r + i]	+= partial_products[ th*nmode*r +  (nmode-2) * r + i] * t->val[it];
				/*
				if(i==0)
					printf("%lf %lf %lf %d\n",  vals[t->ind[nmode-1][it]*r + i], partial_products[(nmode-2) * r + i], t->val[it], t->ind[nmode-1][it]);
				*/
			}
			it++;
	
		}
	}
	mats[nmode-1]->val = vals;

	return 0;
}


int mttkrp_atomic(csf* t, int mode, int r, matrix** mats)
{
	printf("here\n");
	printf("%d %d \n", mode, t->nmode);
	if (mode == (t->nmode)-1)
	{
		return mttkrp_atomic_last(t,mode,r,mats);
	}

	printf("%d %d \n", mode, t->nmode);

	return 0;
}

#endif