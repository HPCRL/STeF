#ifndef MTTKRP_CPP
#define MTTKRP_CPP
#include "../inc/mttkrp.h"



int mttkrp_atomic3(csf* t, int mode, int r, matrix** mats)
{


	TYPE* partial_products;
	int i,j,k,x,y,ind,nmode;
	TYPE* vals;

	nmode = t->nmode;
	partial_products = (TYPE* ) malloc(nmode*r*sizeof(TYPE));
	vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	
	matrix m;


	for(i=0 ; i<mats[mode]->dim1*mats[mode]->dim2 ; i++)
		vals[i] = 0;

	LIKWID_MARKER_INIT;
	LIKWID_MARKER_THREADINIT;
	LIKWID_MARKER_START("Compute");

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
					//if (y == 0)						printf("%lf %lf %d %d\n", partial_products[r + y] * t->val[k] , vals[kk*mats[mode]->dim2 + y], kk*mats[mode]->dim2 + y, kk);
					//printf("%lf\n", t->val[k]);
				}
			}
		}
		
	}

	mats[nmode-1]->val = vals;
	LIKWID_MARKER_STOP("Compute");
	LIKWID_MARKER_CLOSE;
	return 0;
}


int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats, int vec)
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
		th = 0;
	#endif


	printf("num ths %d\n", num_th);
	partial_products = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	inds = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	nnz = t->fiber_count[nmode-1];
	vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	for(i=0 ; i<mats[mode]->dim1*mats[mode]->dim2 ; i++)
		vals[i] = 0;
	

	//printf("nmode is %d nnz is %d \n",nmode,nnz);

	LIKWID_MARKER_INIT;
	LIKWID_MARKER_THREADINIT;
	LIKWID_MARKER_START("Compute");

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

	if(vec == 0)
	{
		for(i=0 ; i<r ; i++)
		{
			vals[t->ind[nmode-1][0]*r + i]	= partial_products[(nmode-2) * r + i] * t->val[0];
		}
	}
	else
	{
		for(i=0 ; i<r ; i++)
		{
			vals[t->ind[nmode-1][0]*r + i]	= partial_products[(nmode-2) * r + i] * t->val[i];
		}	
	}

	for(i = 0 ; i < nmode ; i++)
	{
		inds[i] = 0;
	}

	it = 1;
	//printf("%d\n",nnz);
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
						inds[th*nmode + i] ++;
	
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
						partial_products[th*nmode*r +  ii*r + i] = MAT(mats[ii], t->ind[ii][inds[th*nmode + ii]] ,i) * partial_products[th*nmode*r +  (ii-1)*r + i];
					}
					update ++;
				}
			}
			// use partial products to update last mode
			// Assign all access to a variable
			int xx , yy;
			xx = t->ind[nmode-1][it]*r;
			yy = th*nmode*r +  (nmode-2) * r ;

			if(vec == 0)
			{
				TYPE tval = t->val[it];
				#pragma omp simd
				for(i=0 ; i<r ; i++)
				{
					// put a locking step here
					// This should be atomic
					vals[xx + i]	+= partial_products[yy + i] * tval;
				}
				
			}
			else
			{

				TYPE* tval = (t->val) + it*r;
				#pragma omp simd
				for(i=0 ; i<r ; i++)
				{
					// put a locking step here
					// This should be atomic
					vals[xx + i]	+= partial_products[yy + i] * tval[i];
				}
				//printf("here %lf %lf %lf \n",vals[xx],partial_products[yy],tval[0]);
			}

			it++;
		}
	}
	rem(mats[nmode-1]->val);
	mats[nmode-1]->val = vals;

	LIKWID_MARKER_STOP("Compute");
	LIKWID_MARKER_CLOSE;

	rem(inds);
	rem(partial_products);
	return 0;
}

int mttkrp_fused_init(csf* t,int r)
{
	int i,j;
	long long total_space = 0;
	char* space_sign;

	if(t->intval == NULL)
	{
		t->intval = (TYPE**) malloc((t->nmode)*sizeof(TYPE));
		for(i = 1; i < (t->nmode)-1 ; i++)
		{
			t->intval[i] = (TYPE*) malloc((t->fiber_count[i])*r*sizeof(TYPE));
			if(t->intval[i] == NULL)
			{
				printf("Allocation error in mttkrp fused init\n");
				exit(1);
			}
			else
			{
				total_space += (t->fiber_count[i])*r*sizeof(TYPE);
			}
		}
		if(total_space >= 1073741824)
		{
			// GB
			space_sign = "GB";
			total_space /= 1073741824;
		}
		else if(total_space >= 1048576)
		{
			// MB
			space_sign = "MB";
			total_space /= 1048576;
		}
		else if(total_space >= 1024)
		{
			// KB
			space_sign = "KB";
			total_space /= 1024;
		}
		else
		{
			// B
			space_sign = "B";
		}
		printf("Additional space requirement for the intermediate tensors is %ld%s \n",total_space,space_sign);
	}
	for(i = 1; i < (t->nmode)-1 ; i++)
	{
		for(j=0; j<(t->fiber_count[i])*r ; j++)
		{
			t->intval[i][j] = 0;
		}
	}

	return 0;
}

int mttkrp_atomic_first(csf* t, int mode, int r, matrix** mats)
{
	TYPE *partial_products, *vals;
	int i,ii,it,nmode,nnz;
	idx_t* inds;
	int num_th, th;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		th = 0;
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

	
	LIKWID_MARKER_INIT;
	LIKWID_MARKER_THREADINIT;
	LIKWID_MARKER_START("Compute");


	for(ii = 0; ii < nmode - 1 ; ii++)
	{
		for(i = 0 ; i<r ; i++)
		{
			partial_products[ii * r + i]  =  0 ; //partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][0],i);
		}	
	}

	for(i=0 ; i<r ; i++)
	{
		partial_products[(nmode-2) * r + i] = t->val[0] * MAT(mats[nmode-1],t->ind[nmode-1][0],i);
		//printf("first %lf %lf %lf\n",partial_products[(nmode-2) * r + i], );
	}

	for(i = 0 ; i < nmode ; i++)
	{
		inds[i] = 0;
	}

	it = 1;
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		#ifdef OMP
			th = omp_get_thread_num();
		#endif
		while(it<nnz)
		{
			/*
			printf("partial products ");
			for(i = 0 ; i<(nmode-1)*r; i++)
			{
				printf("%lf ", partial_products[i]);
			}
			printf("\n");
			*/
			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				//inds[th*nmode + nmode-2] ++;
				for(i = nmode-3; i>=0 ; i--)
				{
					if(inds[th*nmode + i+1] + 1 < t->ptr[i][inds[th*nmode + i]+1] )
					{
						update = i+1;
						break;
					}
					/*
					else
					{
						inds[th*nmode + i] ++;
					}
					*/
				}
				/*
				for(i = 0 ; i < nmode ; i++ )
				{
					printf("%d ",inds[i]);
				}
				printf("\n");
				*/
				
				//printf("update %d\n", update);

				for(ii = nmode-2; ii>=update && ii>0 ; ii--)
				{
					int xx = th*nmode*r +  (ii-1)*r;
					int yy = th*nmode*r +  ii*r ;
					int zz = t->ind[ii][inds[th*nmode + ii]];
					#pragma omp simd
					for(i = 0 ; i<r ; i++)
					{
						partial_products[xx + i]  +=   partial_products[yy + i] * MAT(mats[ii],  zz ,i);
					}
				}			
				
				if(update == 0)
				{

					for(i = 0 ; i<r ; i++)
					{
						vals[ (t->ind[0][inds[th*nmode]])*r +i] = partial_products[th*nmode*r + i];
						
						partial_products[th*nmode*r + i] = 0;
						
					}
					update ++;
					inds[ th*nmode ]++;

				}

				for(ii = update; ii <= nmode - 2; ii++)
				{
					for(i = 0 ; i<r ; i++)
					{
						t->intval[ii][inds[th*nmode + ii]*r + i]  =   partial_products[th*nmode*r +  ii*r + i];
						partial_products[th*nmode*r +  ii*r + i] = 0;
					}
					inds[th*nmode + ii] ++;
				}
				/*
				for(i = 0 ; i < nmode ; i++ )
				{
					printf("%d ",inds[i]);
				}
				printf("are inds at exit \n");
				*/

			}
			// use partial products to update last mode
			// Assign all access to a variable
			int xx , yy;
			xx = t->ind[nmode-1][it]*r;
			yy = th*nmode*r +  (nmode-2) * r ;
			TYPE tval = t->val[it];
			TYPE* matval = (mats[nmode-1]->val)+(mats[nmode-1]->dim2)*t->ind[nmode-1][it];
			

			#pragma omp simd
			for(i=0 ; i<r ; i++)
			{
				// put a locking step here
				// This should be atomic
				//vals[xx + i]	+= partial_products[yy + i] * tval[it];
				partial_products[yy + i] += tval * matval[i];
				/*
				if(i==0)
					printf("%lf %lf %lf %d\n",  vals[t->ind[nmode-1][it]*r + i], partial_products[(nmode-2) * r + i], t->val[it], t->ind[nmode-1][it]);
				*/
			}
			it++;
	
		}
	}
	for(ii = nmode-2; ii>0 ; ii--)
	{
		int xx = th*nmode*r +  (ii-1)*r;
		int yy = th*nmode*r +  ii*r ;
		int zz = t->ind[ii][inds[th*nmode + ii]];
		#pragma omp simd
		for(i = 0 ; i<r ; i++)
		{
			partial_products[xx + i]  +=   partial_products[yy + i] * MAT(mats[ii],  zz ,i);
		}
	}			

	{

		for(i = 0 ; i<r ; i++)
		{
			vals[ (t->ind[0][inds[th*nmode]])*r +i] = partial_products[th*nmode*r + i];
			partial_products[th*nmode*r + i] = 0;
			
		}
		inds[ th*nmode ]++;

	}

	for(ii = 1; ii <= nmode - 2; ii++)
	{
		for(i = 0 ; i<r ; i++)
		{
			t->intval[ii][inds[th*nmode + ii]*r + i]  =   partial_products[th*nmode*r +  ii*r + i];
			partial_products[th*nmode*r +  ii*r + i] = 0;
		}
		inds[th*nmode + ii] ++;
	}
	rem(mats[mode]->val);
	mats[mode]->val = vals;


	LIKWID_MARKER_STOP("Compute");
	LIKWID_MARKER_CLOSE;

	rem(inds);
	rem(partial_products);


	return 0;

}

int mttkrp_atomic_middle(csf* t, int mode, int r, matrix** mats)
{
	csf tt = *t;
	tt.nmode = mode+1;
	tt.val = t->intval[mode];
	//tt.nnz = t->fiber_count[mode];
	//printf("%d ",mode);
	/*
	printf("intermediate values\n");
	for(int i=0; i<t->fiber_count[1] ; i++ )
	{
		for(int j=0 ;j<r ; j++)
		{
			printf("%lf ", t->intval[mode][i*r+j]);
		}
		printf("\n");
	}
	*/
	return mttkrp_atomic_last(&tt,mode,r,mats,1);
}

int mttkrp_atomic(csf* t, int mode, int r, matrix** mats)
{
	//printf("here\n");
	//printf("%d %d \n", mode, t->nmode);
	if (mode == 0)
	{
		// return mttkrp_atomic_first(t,mode,mats);

		return mttkrp_atomic_first(t,mode,r,mats);

	}
	else if (mode == (t->nmode)-1)
	{
		return mttkrp_atomic_last(t,mode,r,mats);
	}
	else
	{
		//printf("To be implemented\n");
		return mttkrp_atomic_middle(t,mode,r,mats);
	}


	printf("%d %d \n", mode, t->nmode);

	return 0;
}



#endif