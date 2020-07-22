{
	TYPE* __restrict__ xx , * __restrict__ yy;
	const idx_t row_id = t->ind[nmode-1][it];
	xx = vals + t->ind[nmode-1][it]*r;
	yy = partial_products + (nmode-2) * r ;

	if(vec == 0)
	{
		TYPE tval = t->val[it];
		#ifdef OMP
		mutex_set_lock(mutex,row_id);
		#endif

		#pragma omp simd
		for(int i=0 ; i<r ; i++)
		{
			// put a locking step here
			// This should be atomic
			TYPE increment = yy [i] * tval;
			//printf("before last mode %lf %lf %lf to pos mat[%d][%d][%d] \n",  xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
			//#pragma omp atomic update
			xx [i]	+= increment;
			//printf("last mode %lf %lf %lf to pos mat[%d][%d][%d] \n", xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
		}
		#ifdef OMP
		mutex_unset_lock(mutex,row_id);
		#endif
	}
	else
	{

		TYPE* tval = (t->val) + it*r;
		#ifdef OMP
		mutex_set_lock(mutex,row_id);
		#endif

		#pragma omp simd
		for(int i=0 ; i<r ; i++)
		{
			// put a locking step here
			// This should be atomic
			TYPE increment = yy[i] * tval[i];
			//#pragma omp atomic update
			xx [i]	+= increment;

		}
		#ifdef OMP
		mutex_unset_lock(mutex,row_id);
		#endif
		//printf("here %lf %lf %lf \n",vals[xx],partial_products[yy],tval[0]);
	}
	

	
}