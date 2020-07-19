{
	TYPE * __restrict__ yy;
	//xx = t->ind[nmode-1][it]*r;
	yy = partial_products + (nmode-2) * r ;
	TYPE tval = t->val[it];
	TYPE const * const __restrict__ matval = (mats[nmode-1]->val)+(mats[nmode-1]->dim2)*(t->ind[nmode-1][it]);


	#pragma omp simd
	for(int i=0 ; i<r ; i++)
	{
		// put a locking step here
		// This should be atomic
		//vals[xx + i]	+= partial_products[yy + i] * tval[it];
		yy [i] += tval * matval[i];
		//#pragma omp critical
		//printf("pp %lf tval  %lf matval %lf nnz_id %d and th %d\n",yy[i],tval,matval[i],it,th);
		
	}
	it++;
}