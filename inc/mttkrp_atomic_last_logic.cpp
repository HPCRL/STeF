{
	TYPE* __restrict__ xx , * __restrict__ yy;
	const idx_t row_id = t->ind[nmode-1][it];
	xx = vals + row_id*r;
	yy = partial_products + (nmode-2) * r ;

	//if(vec == 0)
	#ifndef VEC
	{
		//if(it == 1)
		//	printf("saxpy code\n");

		TYPE tval = t->val[it];
		#ifdef OMP
		mutex_set_lock(mutex,row_id);
		#endif
		//if(row_id < 5)		printf("row %d is alloced in %d  %lf\n", row_id, th,xx[0]);
		#pragma omp simd
		for(int i=0 ; i<r ; i++)
		{
			// put a locking step here
			// This should be atomic
			TYPE increment = yy [i] * tval;
			//printf("Increment is %lf\n",Increment );
			//printf("before last mode %lf %lf %lf to pos mat[%d][%d][%d] \n",  xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
			//#pragma omp atomic update
			xx [i]	+= increment;
			//printf("last mode %lf %lf %lf to pos mat[%d][%d][%d] \n", xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
		}
		//if(row_id < 5)		printf("row %d is released in %d val[0] is %lf \n", row_id, th,xx[0]);

		#ifdef OMP
		mutex_unset_lock(mutex,row_id);
		#endif
	}
	#else
	{
		//if(it == 1)
		//	printf("VEC code\n");

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
	#endif

	//printf("first val is %lf %lf \n",vals[0],mats[mode]->val[0]);	
}