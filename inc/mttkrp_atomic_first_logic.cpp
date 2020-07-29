/*
#pragma omp critical
{
for(i = 0 ; i<nmode*r ; i++)
{
	printf("%lf ", partial_products[i]);
}
printf("before updating partial_products th %d\n",th);	
}
*/


for(int ii = nmode-2; ii>=update && ii>0 ; ii--)
{
	TYPE* __restrict__ xx =  partial_products + (ii-1)*r;
	TYPE const * const __restrict__ yy =  partial_products + ii*r ;
	TYPE const * const __restrict__ zz = (mats[ii]->val) + (t->ind[ii][inds[ii]])*(mats[ii]->dim2);
	#pragma omp simd
	for(int i = 0 ; i<r ; i++)
	{
		// partial_products[xx + i]  +=   partial_products[yy + i] * MAT(mats[ii],  zz ,i);
		xx[i] += yy[i] * zz[i];
	}
}		
/*
#pragma omp critical
{
for(i = 0 ; i<nmode*r ; i++)
{
	printf("%lf ", partial_products[i]);
}
printf("after updating partial_products th %d\n",th);	
}
*/

if(update == 0)
{
	if(DOT_PARALLEL_DEPTH <= 1)
	{
		TYPE * __restrict__ xx = vals + (t->ind[0][inds[0]])*r;
		TYPE const * const __restrict__ yy = partial_products;
		#pragma omp simd
		for(int i = 0 ; i<r ; i++)
		{
			xx[i] = yy[i];
		}
	}
	else
	{	
		TYPE* __restrict__ xx = vals + (t->ind[0][inds[0]])*r;
		TYPE const * const __restrict__ yy = partial_products;
		#ifdef OMP
		const idx_t row_id = 0;
		mutex_set_lock(mutex,row_id);
		//printf("here\n");
		#endif

		#pragma omp simd
		for(int i = 0 ; i<r ; i++)
		{
			//#pragma omp atomic update
			xx[i] += yy[i];
		}

		#ifdef OMP 
		mutex_unset_lock(mutex,row_id);
		#endif
	}

	#pragma omp simd
	for(int i = 0 ; i<r ; i++)
	{	
		//printf("partial_result at %d added to output at %d %lf in th %d\n",i,(t->ind[0][inds[0]])*r + i,partial_products[i] ,th);
		partial_products[ i] = 0;
		
	}
	update ++;
	inds[ 0]++;

}

for(int ii = update; ii < nmode - 1; ii++)
{
	if(DOT_PARALLEL_DEPTH > ii+1)
	{
		TYPE* __restrict__ xx = t->intval[ii]  + (inds[ii]*r);

		TYPE* __restrict__ yy = partial_products + ii*r;
		//const idx_t row_id = 

		#ifdef OMP
		const idx_t row_id = ii;
		mutex_set_lock(mutex,row_id);
		#endif

		#pragma omp simd
		for(int i = 0 ; i<r ; i++)
		{
			//#pragma omp atomic update
			xx[i] += yy[i];
		}

		#ifdef OMP 
		mutex_unset_lock(mutex,row_id);
		#endif

	}
	else
	{	
		#pragma omp simd 
		for(int i = 0 ; i<r ; i++)
		{
			t->intval[ii][inds[ii]*r + i]  =   partial_products[ ii*r + i];			
		}
	}	
	
	#pragma omp simd
	for(int i = 0 ; i<r ; i++)
	{
		//printf("partial_result at %d added to intval at %d %lf in th %d\n",ii*r+i,inds[ii]*r+i,partial_products[ii*r+i] ,th);
		partial_products[ ii*r + i] = 0;		
	}		
	inds[ ii] ++;
}