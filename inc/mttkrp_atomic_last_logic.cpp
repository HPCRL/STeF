{
	TYPE* xx , *yy;
	xx = vals + t->ind[nmode-1][it]*r;
	yy = partial_products + (nmode-2) * r ;

	if(vec == 0)
	{
		TYPE tval = t->val[it];
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
		
	}
	else
	{

		TYPE* tval = (t->val) + it*r;
		#pragma omp simd
		for(int i=0 ; i<r ; i++)
		{
			// put a locking step here
			// This should be atomic
			TYPE increment = yy[i] * tval[i];
			//#pragma omp atomic update
			xx [i]	+= increment;

		}
		//printf("here %lf %lf %lf \n",vals[xx],partial_products[yy],tval[0]);
	}
	

	
}