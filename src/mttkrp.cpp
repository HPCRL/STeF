#ifndef MTTKRP_CPP
#define MTTKRP_CPP
#include "../inc/mttkrp.h"

#ifdef OMP
mutex_array* mutex = NULL;
#endif

int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats, int vec, int profile)
{
	/* 
	if(t->nmode == 3)
	{
		return mttkrp_atomic3(t,mode,r,mats);
	}
	*/



	TYPE* partial_products_all;
	int nmode;
	idx_t* inds_all;
	int num_th;
	TYPE* temp_res_all;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		int th = 0;
	#endif
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	temp_res_all = (TYPE* ) malloc(num_th*(r+64)*sizeof(TYPE));
	inds_all = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	long long int nnz = t->fiber_count[nmode-1];
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	TYPE* vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//printf("nmode is %d nnz is %d \n",nmode,nnz);
	
	if(profile == mode)
	{
		printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}
	

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if (profile == mode)
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}

	// it = 0;
	// DO the process for the first nnz
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
			int th = omp_get_thread_num();
			if(VERBOSE == VERBOSE_DEBUG)
			printf("th id is %d\n",th);
		#endif

		idx_t  end;
		long long int it = (th*nnz)/num_th;
		end = ((th+1)*nnz)/num_th;
		if(VERBOSE == VERBOSE_DEBUG)
			printf("start %lld end is %d for thread %d \n",it,end,th);

		TYPE* partial_products;	
		idx_t* inds = inds_all + th*nmode;
		partial_products = partial_products_all + th*nmode*r;
		//TYPE* temp_res = temp_res_all + th*(r+64);
		//TYPE* vals;
		/*
		if(num_th == 1)
		{
			//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
			vals = mats[mode]->val;
			
			
		}
		else // Use private copies
		{	
			vals = t->private_mats[th]->val;
			//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		}
		memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		*/

		auto time_start = std::chrono::high_resolution_clock::now();
		find_inds(inds,t,it);
		
		
	
		for(int i = 0 ; i<r ; i++)
		{
			partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
		}
	
		for(int ii = 1; ii < nmode - 1 ; ii++)
		{
			for(int i = 0 ; i<r ; i++)
			{
				partial_products[ii * r + i]  = partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][inds[ii]],i);
			}	
		}
	
		#include "../inc/mttkrp_atomic_last_logic.cpp"
		it ++;
		while(it<end)
		{

			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				inds[nmode-2] ++;
				for(int i = nmode-3; i>=0 ; i--)
				{
					if(inds[i+1]  < t->ptr[i][inds[i]+1] )
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
					for(int i = 0 ; i<r ; i++)
					{
						partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
						//printf("update %lf %lf \n, ",MAT(mats[0], t->ind[0][inds[nmode*th  +0]] ,i), partial_products[th*nmode*r + i]);
					}
					update ++;
				}
	
	
	
				for(int ii = update; ii<nmode-1; ii++)
				{
					TYPE*  xx = partial_products + ii*r;
					TYPE  *  yy = (mats[ii]->val) + (t->ind[ii][inds[ii]])*(mats[ii]->dim2);
					TYPE  *  zz = partial_products + (ii-1)*r;
					for(int i = 0 ; i<r ; i++)
					{
						//partial_products[ ii*r + i] = MAT(mats[ii], t->ind[ii][inds[ii]] ,i) * partial_products[ (ii-1)*r + i];
						xx[i] = yy[i] * zz[i];
					}
					update ++;
				}
			}
			// use partial products to update last mode
			// Assign all access to a variable
			//TYPE* xx , *yy;
			//xx = vals + t->ind[nmode-1][it]*r;
			//yy = partial_products + (nmode-2) * r ;

			#include "../inc/mttkrp_atomic_last_logic.cpp"

			it++;
		}
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Atomic kernel time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
	}

	//rem(mats[nmode-1]->val);
	//mats[nmode-1]->val = vals;
	//printf("here first val is %lf %lf \n",vals[0],mats[mode]->val[0]);

	LIKWID_MARKER_CLOSE;
	t->num_th = num_th;
	/*
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		//reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Atomic kernel time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}
	*/
	rem(inds_all);
	rem(partial_products_all);
	rem(temp_res_all);
	return 0;
}

int mttkrp_private_last(csf* t, int mode, int r, matrix** mats, int vec, int profile)
{
	TYPE* partial_products_all;
	int nmode;
	idx_t* inds_all;
	int num_th;
	TYPE* temp_res_all;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		int th = 0;
	#endif
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	temp_res_all = (TYPE* ) malloc(num_th*(r+64)*sizeof(TYPE));
	inds_all = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	long long int nnz = t->fiber_count[nmode-1];
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//TYPE* vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//printf("nmode is %d nnz is %d \n",nmode,nnz);
	
	if(profile == mode)
	{
		printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}
	

	// Set matrices to 0 here
	if(num_th == 1)
	{
		memset(mats[mode]->val, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	}
	else
	{
		#pragma omp parallel
		for(int i =0 ; i<num_th ; i++)
		{
			TYPE* val0 = t->private_mats[i]->val;
			long long int range = ((long long int) mats[mode]->dim1)*(mats[mode]->dim2);
			int steps = num_th*4;
			long long int jump = range/steps;
			#pragma omp for schedule (dynamic,1)
			for(int j=0 ; j < steps ; j++)
			{
				memset(val0 + j*jump , 0 , jump*sizeof(TYPE));
			}
			#pragma omp single
			memset(val0 + steps*jump, 0 , (range - steps*jump)*sizeof(TYPE));
		}
		
			
	}

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if (profile == mode)
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}

	// it = 0;
	// DO the process for the first nnz
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
			int th = omp_get_thread_num();
			if(VERBOSE == VERBOSE_DEBUG)
			printf("th id is %d\n",th);
		#endif

		idx_t  end;
		long long int it = (th*nnz)/num_th;
		end = ((th+1)*nnz)/num_th;
		if(VERBOSE == VERBOSE_DEBUG)
			printf("start %lld end is %d for thread %d \n",it,end,th);

		TYPE* partial_products;	
		idx_t* inds = inds_all + th*nmode;
		partial_products = partial_products_all + th*nmode*r;
		TYPE* temp_res = temp_res_all + th*(r+64);
		TYPE* vals;
		
		if(num_th == 1)
		{
			//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
			vals = mats[mode]->val;
			
			
		}
		else // Use private copies
		{	
			vals = t->private_mats[th]->val;
			//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		}
		//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		

		auto time_start = std::chrono::high_resolution_clock::now();
		find_inds(inds,t,it);
		
		
	
		for(int i = 0 ; i<r ; i++)
		{
			partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
		}
	
		for(int ii = 1; ii < nmode - 1 ; ii++)
		{
			for(int i = 0 ; i<r ; i++)
			{
				partial_products[ii * r + i]  = partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][inds[ii]],i);
			}	
		}
	
		#include "../inc/mttkrp_atomic_last_logic.cpp"
		it ++;
		while(it<end)
		{

			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				inds[nmode-2] ++;
				for(int i = nmode-3; i>=0 ; i--)
				{
					if(inds[i+1]  < t->ptr[i][inds[i]+1] )
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
					for(int i = 0 ; i<r ; i++)
					{
						partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
						//printf("update %lf %lf \n, ",MAT(mats[0], t->ind[0][inds[nmode*th  +0]] ,i), partial_products[th*nmode*r + i]);
					}
					update ++;
				}
	
	
	
				for(int ii = update; ii<nmode-1; ii++)
				{
					TYPE*  xx = partial_products + ii*r;
					TYPE  *  yy = (mats[ii]->val) + (t->ind[ii][inds[ii]])*(mats[ii]->dim2);
					TYPE  *  zz = partial_products + (ii-1)*r;
					for(int i = 0 ; i<r ; i++)
					{
						//partial_products[ ii*r + i] = MAT(mats[ii], t->ind[ii][inds[ii]] ,i) * partial_products[ (ii-1)*r + i];
						xx[i] = yy[i] * zz[i];
					}
					update ++;
				}
			}
			// use partial products to update last mode
			// Assign all access to a variable
			//TYPE* xx , *yy;
			//xx = vals + t->ind[nmode-1][it]*r;
			//yy = partial_products + (nmode-2) * r ;

			#include "../inc/mttkrp_atomic_last_logic.cpp"

			it++;
		}
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Private kernel time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
	}

	//rem(mats[nmode-1]->val);
	//mats[nmode-1]->val = vals;
	//printf("here first val is %lf %lf \n",vals[0],mats[mode]->val[0]);

	LIKWID_MARKER_CLOSE;
	t->num_th = num_th;
	
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Private kernel time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}
	
	rem(inds_all);
	rem(partial_products_all);
	rem(temp_res_all);
	return 0;
}

int mttkrp_atomic_last_vec(csf* t, int mode, int r, matrix** mats, int vec, int profile)
{
	/* 
	if(t->nmode == 3)
	{
		return mttkrp_atomic3(t,mode,r,mats);
	}
	*/
	#define VEC


	TYPE* partial_products_all;
	int nmode;
	idx_t* inds_all;
	int num_th;
	TYPE* temp_res_all;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		int th = 0;
	#endif
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	temp_res_all = (TYPE* ) malloc(num_th*(r+64)*sizeof(TYPE));
	inds_all = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	long long int nnz = t->fiber_count[nmode-1];
	//TYPE* vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	TYPE* vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//printf("nmode is %d nnz is %d \n",nmode,nnz);
	
	if(profile == mode)
	{
		printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}
	

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if (profile == mode)
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}

	// it = 0;
	// DO the process for the first nnz
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
			int th = omp_get_thread_num();
			if(VERBOSE == VERBOSE_DEBUG)
			printf("th id is %d\n",th);
		#endif

		idx_t it, end;
		it = (th*nnz)/num_th;
		end = ((th+1)*nnz)/num_th;
		if(VERBOSE == VERBOSE_DEBUG)
			printf("start %d end is %d for thread %d \n",it,end,th);

		TYPE* partial_products;	
		idx_t* inds = inds_all + th*nmode;
		partial_products = partial_products_all + th*nmode*r;
		//TYPE* temp_res = temp_res_all + th*(r+64);
		/*
		TYPE* vals;

		if(num_th == 1)
		{
			//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
			vals = mats[mode]->val;
			
			
		}
		else // Use private copies
		{	
			vals = t->private_mats[th]->val;
			//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		}
		memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		*/
		auto time_start = std::chrono::high_resolution_clock::now();
		find_inds(inds,t,it);
		
		
	
		for(int i = 0 ; i<r ; i++)
		{
			partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
		}
	
		for(int ii = 1; ii < nmode - 1 ; ii++)
		{
			for(int i = 0 ; i<r ; i++)
			{
				partial_products[ii * r + i]  = partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][inds[ii]],i);
			}	
		}
	
		#include "../inc/mttkrp_atomic_last_logic.cpp"
		it ++;
		while(it<end)
		{

			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				inds[nmode-2] ++;
				for(int i = nmode-3; i>=0 ; i--)
				{
					if(inds[i+1]  < t->ptr[i][inds[i]+1] )
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
					for(int i = 0 ; i<r ; i++)
					{
						partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
						//printf("update %lf %lf \n, ",MAT(mats[0], t->ind[0][inds[nmode*th  +0]] ,i), partial_products[th*nmode*r + i]);
					}
					update ++;
				}
	
	
	
				for(int ii = update; ii<nmode-1; ii++)
				{
					TYPE*  xx = partial_products + ii*r;
					TYPE  *  yy = (mats[ii]->val) + (t->ind[ii][inds[ii]])*(mats[ii]->dim2);
					TYPE  *  zz = partial_products + (ii-1)*r;
					for(int i = 0 ; i<r ; i++)
					{
						//partial_products[ ii*r + i] = MAT(mats[ii], t->ind[ii][inds[ii]] ,i) * partial_products[ (ii-1)*r + i];
						xx[i] = yy[i] * zz[i];
					}
					update ++;
				}
			}
			// use partial products to update last mode
			// Assign all access to a variable
			//TYPE* xx , *yy;
			//xx = vals + t->ind[nmode-1][it]*r;
			//yy = partial_products + (nmode-2) * r ;

			#include "../inc/mttkrp_atomic_last_logic.cpp"

			it++;
		}
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Atomic kernel time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
	}

	//rem(mats[nmode-1]->val);
	//mats[nmode-1]->val = vals;

	LIKWID_MARKER_CLOSE;
	t->num_th = num_th;
	/*if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Atomic kernel time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}*/

	rem(inds_all);
	rem(partial_products_all);
	rem(temp_res_all);
	#undef VEC
	return 0;
}


int mttkrp_private_last_vec(csf* t, int mode, int r, matrix** mats, int vec, int profile)
{
	#define VEC
	TYPE* partial_products_all;
	int nmode;
	idx_t* inds_all;
	int num_th;
	TYPE* temp_res_all;

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		int th = 0;
	#endif
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	temp_res_all = (TYPE* ) malloc(num_th*(r+64)*sizeof(TYPE));
	inds_all = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	long long int nnz = t->fiber_count[nmode-1];
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//TYPE* vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//printf("nmode is %d nnz is %d \n",nmode,nnz);
	
	if(profile == mode)
	{
		printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}
	

	// Set matrices to 0 here
	if(num_th == 1)
	{
		memset(mats[mode]->val, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	}
	else
	{
		#pragma omp parallel
		for(int i =0 ; i<num_th ; i++)
		{
			TYPE* val0 = t->private_mats[i]->val;
			long long int range = ((long long int) mats[mode]->dim1)*(mats[mode]->dim2);
			int steps = num_th*4;
			long long int jump = range/steps;
			#pragma omp for schedule (dynamic,1)
			for(int j=0 ; j < steps ; j++)
			{
				memset(val0 + j*jump , 0 , jump*sizeof(TYPE));
			}
			#pragma omp single
			memset(val0 + steps*jump, 0 , (range - steps*jump)*sizeof(TYPE));
		}
		
	}

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if (profile == mode)
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}

	// it = 0;
	// DO the process for the first nnz
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
			int th = omp_get_thread_num();
			if(VERBOSE == VERBOSE_DEBUG)
			printf("th id is %d\n",th);
		#endif

		idx_t  end;
		long long int it = (th*nnz)/num_th;
		end = ((th+1)*nnz)/num_th;
		if(VERBOSE == VERBOSE_DEBUG)
			printf("start %lld end is %d for thread %d \n",it,end,th);

		TYPE* partial_products;	
		idx_t* inds = inds_all + th*nmode;
		partial_products = partial_products_all + th*nmode*r;
		TYPE* temp_res = temp_res_all + th*(r+64);
		TYPE* vals;
		
		if(num_th == 1)
		{
			//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
			vals = mats[mode]->val;
			
			
		}
		else // Use private copies
		{	
			vals = t->private_mats[th]->val;
			//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		}
		//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
		

		auto time_start = std::chrono::high_resolution_clock::now();
		find_inds(inds,t,it);
		
		
	
		for(int i = 0 ; i<r ; i++)
		{
			partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
		}
	
		for(int ii = 1; ii < nmode - 1 ; ii++)
		{
			for(int i = 0 ; i<r ; i++)
			{
				partial_products[ii * r + i]  = partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][inds[ii]],i);
			}	
		}
	
		#include "../inc/mttkrp_atomic_last_logic.cpp"
		it ++;
		while(it<end)
		{

			// traverse nnz
			// if there is a new fiber
			// --update partial_products
			
			if(it == t->ptr[nmode-2][inds[nmode-2]+1])
			{
				int update = 0;
				inds[nmode-2] ++;
				for(int i = nmode-3; i>=0 ; i--)
				{
					if(inds[i+1]  < t->ptr[i][inds[i]+1] )
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
					for(int i = 0 ; i<r ; i++)
					{
						partial_products[i] = MAT(mats[0], t->ind[0][inds[0]] ,i);
						//printf("update %lf %lf \n, ",MAT(mats[0], t->ind[0][inds[nmode*th  +0]] ,i), partial_products[th*nmode*r + i]);
					}
					update ++;
				}
	
	
	
				for(int ii = update; ii<nmode-1; ii++)
				{
					TYPE*  xx = partial_products + ii*r;
					TYPE  *  yy = (mats[ii]->val) + (t->ind[ii][inds[ii]])*(mats[ii]->dim2);
					TYPE  *  zz = partial_products + (ii-1)*r;
					for(int i = 0 ; i<r ; i++)
					{
						//partial_products[ ii*r + i] = MAT(mats[ii], t->ind[ii][inds[ii]] ,i) * partial_products[ (ii-1)*r + i];
						xx[i] = yy[i] * zz[i];
					}
					update ++;
				}
			}
			// use partial products to update last mode
			// Assign all access to a variable
			//TYPE* xx , *yy;
			//xx = vals + t->ind[nmode-1][it]*r;
			//yy = partial_products + (nmode-2) * r ;

			#include "../inc/mttkrp_atomic_last_logic.cpp"

			it++;
		}
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Private kernel time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
	}

	//rem(mats[nmode-1]->val);
	//mats[nmode-1]->val = vals;
	//printf("here first val is %lf %lf \n",vals[0],mats[mode]->val[0]);

	LIKWID_MARKER_CLOSE;
	t->num_th = num_th;
	
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Private kernel time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}
	
	rem(inds_all);
	rem(partial_products_all);
	rem(temp_res_all);
	#undef VEC
	return 0;
}


int mttkrp_fused_init(csf* t,int r)
{
	int i,j;
	long long total_space = 0;
	char* space_sign;

	#ifdef OMP
	if (mutex == NULL)
	{
		mutex = mutex_alloc_custom((t->mlen)[t->nmode-1] , 16);
		//mutex = mutex_alloc_custom(1024 , 16); // This is what splatt is using
	}

	int num_th = omp_get_max_threads();
	t->private_mats = (matrix** ) malloc(num_th*sizeof(matrix*));

	idx_t max_len = t->mlen[0];
	for(int i=1; i<t->nmode ; i++)
	{	
		if(t->mlen[i] > max_len)
			max_len = t->mlen[i];
	}
	for(int i=0; i<num_th ; i++)
	{
		t->private_mats[i] = create_matrix(max_len, r, 0);
	}
	t->num_th = num_th;
	total_space += max_len*r*num_th;
	#else
	t->num_th = 1;
	#endif


	if(t->intval == NULL)
	{
		t->intval = (TYPE**) malloc((t->nmode)*sizeof(TYPE));
		t->intval[0] = NULL;
		t->intval[(t->nmode)-1] = NULL;
		for(i = 1; i < (t->nmode)-1 ; i++)
		{
			t->intval[i] = (TYPE*) malloc((t->fiber_count[i])*r*sizeof(TYPE));
			if(t->intval[i] == NULL)
			{
				printf("SpTL ERROR: Allocation error in mttkrp fused init\n");
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
		printf("Additional space requirement for the intermediate tensors is %llu%s \n",total_space,space_sign);
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


int mttkrp_fused_free(csf* t,int r)
{
	return 0;
}
// Finds the first nnz id corresponding to fiber pointed by an index
idx_t find_nnz_pos(csf* t, int depth, idx_t loc) 
{
	int nmode = t->nmode;
	idx_t nnz_pos = loc;
	int mode_id = depth-1;
	while ( mode_id < nmode-1)
	{
		nnz_pos = t->ptr[mode_id][nnz_pos];
		mode_id ++;
	}
	
	return nnz_pos;	
}

int dist_dot_work(idx_t* inds ,csf* t,int p,idx_t* count, int th,int depth)
{
	int nmode = t->nmode;
	long long int nnz = t->fiber_count[nmode-1];
	//idx_t loc = 0;
	idx_t start = 0;
	idx_t end = 0;
	idx_t goal = (th*nnz)/p;
	idx_t end_goal = ((th+1)*nnz)/p;

	if(th == 0)
	{
		start = 0;
	}
	else
	{
		int bstart = 0;
		int bend = t->fiber_count[depth-1];

		// DO a binary search to find the location

		while (bstart < bend - 1)
		{
			int pivot = (bstart+bend)/2;
			int pos = find_nnz_pos(t,depth,pivot);
			if (pos < goal)
			{
				bstart = pivot;
			}
			else if (pos > goal)
			{
				bend = pivot;
			}
			else
			{
				bstart = pivot;
				bend = pivot + 1;
			}
		}
		start = bstart;
	}

	
	// Starting position for the thread found

	// Now find the ending position for the thread
	if(end_goal == nnz)
	{
		end = t->fiber_count[depth-1];
	}
	else
	{
		int bstart = 0;
		int bend = t->fiber_count[depth-1];

		// DO a binary search to find the location

		while (bstart < bend - 1)
		{
			int pivot = (bstart+bend)/2;
			int pos = find_nnz_pos(t,depth,pivot);
			if (pos < end_goal)
			{
				bstart = pivot;
			}
			else if (pos > end_goal)
			{
				bend = pivot;
			}
			else
			{
				bstart = pivot;
				bend = pivot + 1;
			}
		}
		end = bstart;
	}
	
	

	if(start >= end)
	{
		printf("SPTL ERROR: No work is available for thread %d\n", th );
		return 1;
	}
	else
	{
		//auto start_time = std::chrono::high_resolution_clock::now();
		idx_t start_pos = find_nnz_pos(t,depth,start);
		
		idx_t end_pos = find_nnz_pos(t,depth,end);
		*count = end_pos-start_pos;
		//printf("start pos is %d\n",start_pos );
		find_inds(inds,t,start_pos);

		//auto end_time = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> diff = end_time-start_time;
		//total += diff.count();
		//printf("preprocessing %lf \n",diff.count());
		return 0;
	}

}

int mttkrp_atomic_first(csf* t, int mode, int r, matrix** mats, int profile)
{
	TYPE *partial_products_all, *vals;
	int nmode;
	idx_t* inds_all;
	int num_th;

	nmode = t->nmode;

	#ifdef OMP
		num_th = omp_get_max_threads();

	#else
		num_th = 1;
		int th = 0;
	#endif

	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*nmode*r*sizeof(TYPE));
	inds_all = (idx_t* ) malloc(num_th * nmode* sizeof(idx_t));
	//int nnz = t->fiber_count[nmode-1];
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	vals = mats[mode]->val;
	//#pragma omp parallel
	for(int i=0 ; i<mats[mode]->dim1*mats[mode]->dim2 ; i++)
		vals[i] = 0;




	// it = 0;
	// DO the process for the first nnz

	if(profile == mode)
	{
		printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}
		
		
	
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if (profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
			int th = omp_get_thread_num();
		#endif


		TYPE* partial_products;	
		idx_t* inds = inds_all + th*nmode;
		partial_products = partial_products_all + th*nmode*r;
		idx_t num_it;

		
		
		dist_dot_work(inds,t,num_th,&num_it,th);
			
		
		if(VERBOSE == VERBOSE_DEBUG)
			printf("first nnz for thread %d is %d and nnz count is %d\n",th, inds[nmode-1], num_it);
		idx_t it_start = inds[nmode-1];
		idx_t it = it_start;

		auto time_start = std::chrono::high_resolution_clock::now();
		for(int ii = 0; ii < nmode - 1 ; ii++)
		{
			#pragma omp simd
			for(int i = 0 ; i<r ; i++)
			{
				partial_products[ii * r + i]  =  0 ; //partial_products[(ii-1) * r + i] * MAT(mats[ii], t->ind[ii][0],i);
			}	
		}
	
		{
			#include "../inc/mttkrp_atomic_first_nnz_logic.cpp"
			//inds[nmode-1] ++ ;
		}
	
		
		{
		
			while(it < num_it+it_start)
			{
				// traverse nnz
				// if there is a new fiber
				// --update partial_products
				
				if(it == t->ptr[nmode-2][inds[nmode-2]+1])
				{
					int update = 0;
					//inds[th*nmode + nmode-2] ++;
					for(int i = nmode-3; i>=0 ; i--)
					{
						if(inds[i+1] + 1 < t->ptr[i][inds[i]+1] )
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
	
					#include "../inc/mttkrp_atomic_first_logic.cpp"
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
				#include "../inc/mttkrp_atomic_first_nnz_logic.cpp"
				//inds[nmode-1] ++ ;
		
			}
		}


	
		int update = 0;
		#include "../inc/mttkrp_atomic_first_logic.cpp"

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Atomic kernel time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());	

		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}

	//rem(mats[mode]->val);
	mats[mode]->val = vals;

	
	LIKWID_MARKER_CLOSE;
	

	rem(inds_all);
	rem(partial_products_all);


	return 0;

}

int mttkrp_atomic_middle(csf* t, int mode, int r, matrix** mats, int profile)
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
	return mttkrp_atomic_last_vec(&tt,mode,r,mats,1,profile);
}

int mttkrp_private_middle(csf* t, int mode, int r, matrix** mats, int profile)
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
	return mttkrp_private_last_vec(&tt,mode,r,mats,1,profile);
}


int mttkrp_atomic(csf* t, int mode, int r, matrix** mats, int profile)
{
	//printf("here\n");
	//printf("%d %d \n", mode, t->nmode);
	if (mode == 0)
	{
		// return mttkrp_atomic_first(t,mode,mats);

		return mttkrp_atomic_first(t,mode,r,mats, profile);
		//return mttkrp_hardwired_first(t,mode,r,mats, profile);

	}
	else if (mode == (t->nmode)-1)
	{
		return mttkrp_atomic_last(t,mode,r,mats,0 , profile);
	}
	else
	{
		//printf("To be implemented\n");
		return mttkrp_atomic_middle(t,mode,r,mats, profile);
	}


	printf("%d %d \n", mode, t->nmode);

	return 0;
}

int mttkrp_private(csf* t, int mode, int r, matrix** mats, int profile)
{
	//printf("here\n");
	//printf("%d %d \n", mode, t->nmode);
	if (mode == 0)
	{
		// return mttkrp_atomic_first(t,mode,mats);

		return mttkrp_atomic_first(t,mode,r,mats, profile);
		//return mttkrp_hardwired_first(t,mode,r,mats, profile);

	}
	else if (mode == (t->nmode)-1)
	{
		return mttkrp_private_last(t,mode,r,mats,0 , profile);
	}
	else
	{
		//printf("To be implemented\n");
		return mttkrp_private_middle(t,mode,r,mats, profile);
	}


	printf("%d %d \n", mode, t->nmode);

	return 0;
}

int mttkrp_test(coo* dt, int mode, int r, matrix** mats)
{
	idx_t i,j,y;
	idx_t nnz = dt -> nnz;
	idx_t size = (mats[mode]->dim1)*(mats[mode]->dim2);
	int nmode = dt -> nmode;
	TYPE* vals = (TYPE* ) malloc(size*sizeof(TYPE));
	TYPE* accum = (TYPE* ) malloc(r*sizeof(TYPE));
	for(i=0; i<size ; i++)
	{
		vals[i] = 0;
	}	

	for(i=0; i<nnz ; i++)
	{
		for(y = 0; y<r ; y++)
			accum[y] = dt -> val[i];  // Load Tensor values 
		/*
		printf("accum is %lf\n", accum[0]);
		for(int m = 0; m < nmode ; m++)
			printf("%d ", dt->ind[i*nmode + m] );
		printf("\n");
		*/

		for(int m = 0; m < nmode ; m++)
		{
			if(m == mode )
			{
				continue;
			}
			else
			{
				idx_t mode_id = m;
				idx_t* ptr = dt->ind + i*nmode;
				idx_t row_id = ptr[ mode_id ];
				//printf("row_id is %d mode_id is %d\n",row_id,mode_id);
				
				TYPE const * const __restrict__ inrow = mats[m]->val + row_id*(mats[m]->dim2);
				#pragma omp simd
				for(y=0 ; y<r ; y++)
				{
					accum[y] *= inrow[y]; // Multiply tensor with the dense matrix
				}
				//printf("accum is %lf\n", accum[0]);
			}

		}
		idx_t mode_id = mode;
		idx_t* ptr = dt->ind + i*nmode;
		idx_t row_id = ptr[ mode_id ];
		TYPE * inrow = vals + row_id*(mats[mode]->dim2);
		for(y=0 ; y<r ; y++)
		{
			inrow[y] += accum[y]; // Multiply tensor with the dense matrix
		}
		//printf("accum is %lf inrow is %lf row is %d\n", accum[0],inrow[0],row_id);
	}

	int num_diff = 0;
	TYPE sum_base= 0, sum_mat = 0;

	for(i = 0 ; i<mats[mode]->dim1; i++)
	{
		for(j = 0 ; j<r ; j++)
		{
			idx_t pos = (i * mats[mode]->dim2) + j;
			TYPE old_val = mats[mode]-> val [pos];
			TYPE new_val = vals[pos];
			sum_base += new_val;
			sum_mat += old_val;

			TYPE diff = old_val - new_val;
			diff /= new_val;
			if(diff < 0 )
				diff = -diff;


			if (diff > CORRECTNESS_THRESHOLD)
			{
				num_diff++;
				if(VERBOSE == VERBOSE_DEBUG || num_diff < 10)
					printf("Correctness Error at mats[%d][%d][%d]. val is %lf, but must be %lf\n", mode,i,j,old_val,new_val );
			}
		}
	}

	printf("total diff is %d / %d. Sums are %lf and %lf \n", num_diff , size, sum_mat, sum_base);
	//rem(vals);
	rem(accum);
	rem(vals);
	return 0;
}


int mttkrp_hardwired(csf* t, int mode, int r, matrix** mats, int profile)
{
	if (mode == 0)
	{
		return mttkrp_hardwired_first(t,mode,r,mats, profile);
		//return mttkrp_hardwired_first(t,mode,r,mats, profile);

	}
	else if (mode == (t->nmode)-1)
	{
		return mttkrp_hardwired_last(t,mode,r,mats, profile);
	}
	else
	{
		//printf("To be implemented\n");
		return mttkrp_hardwired_middle(t,mode,r,mats, profile);
	}


	printf("%d %d \n", mode, t->nmode);

	return 0;
}


int mttkrp_hardwired_first(csf* t, int mode, int r, matrix** mats, int profile)
{

	
	int nmode = t->nmode;
	if (nmode == 3)
		return mttkrp_hardwired_first_3(t,mode,r,mats,profile);
	else if (nmode == 4)
		return mttkrp_hardwired_first_4(t,mode,r,mats,profile);
	else if (nmode == 5)
		return mttkrp_hardwired_first_5(t,mode,r,mats,profile);
	else
		return 1;

	//return 0;

}

int mttkrp_hardwired_last(csf* t, int mode, int r, matrix** mats, int profile)
{
	mutex_array* arr = NULL;
	#ifdef OMP 
	arr = mutex;
	#endif

	int nmode = t->nmode;
	if(nmode == 3)
	{
		return mttkrp_hardwired_last_3(t,mode,r,mats,arr,profile);
	}
	else if(nmode == 4)
	{
		return mttkrp_hardwired_last_4(t,mode,r,mats,arr,profile);
	}
	else if(nmode == 5)
	{
		return mttkrp_hardwired_last_5(t,mode,r,mats,arr,profile);
	}
	return 1;
}


int mttkrp_hardwired_middle(csf* t, int mode, int r, matrix** mats, int profile)
{
	csf tt = *t;
	tt.nmode = mode+1;
	tt.val = t->intval[mode];

	mutex_array* arr = NULL;
	#ifdef OMP 
	arr = mutex;
	#endif

	if(tt.nmode == 2)
	{
		return mttkrp_hardwired_last_vec_2(t,mode,r,mats,arr,profile);
	}
	else if(tt.nmode == 3)
	{
		return mttkrp_hardwired_last_vec_3(t,mode,r,mats,arr,profile);
	}
	else if(tt.nmode == 4)
	{
		return mttkrp_hardwired_last_vec_4(t,mode,r,mats,arr,profile);
	}

	return 1;
}


#endif
