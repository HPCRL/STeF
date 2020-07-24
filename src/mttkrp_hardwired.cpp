#ifndef MTTKRP_HARDWIRED_CPP
#define MTTKRP_HARDWIRED_CPP
#include "../inc/mttkrp_hardwired.h"

int reduce(csf* t, int r, matrix* mat)
{
	memset (mat->val, 0 , sizeof(TYPE)*(mat->dim1)*(mat->dim2));
	int parallel_threshold = 1000;
	for(int i=0; i<t->num_th ; i++)
	{
		if (mat->dim1 > parallel_threshold)
		{
			#pragma omp parallel for
			for(int row_id = 0; row_id < mat->dim1 ; row_id++)
			{
				TYPE* outval = mat->val + row_id * (mat->dim2);
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
		else
		{
			for(int row_id = 0; row_id < mat->dim1 ; row_id++)
			{
				TYPE* outval = mat->val + row_id * (mat->dim2);
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

int mttkrp_hardwired_first_3(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD;
	#ifdef OMP
	num_th = omp_get_max_threads();
	//printf("here \n");
	#endif

	printf("num ths %d\n", num_th);

	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

	for(int i = 0 ; i < num_th*partial_results_size ; i++)
		partial_results_all[i] = 0;

	set_matrix(*mats[0],0);

	#ifdef OMP
	#pragma omp parallel 
	#endif
	{
		int th = 0;
		#ifdef OMP
		th = omp_get_thread_num();
		#endif
		auto time_start = std::chrono::high_resolution_clock::now();
		TYPE* partial_results = partial_results_all + th*partial_results_size;
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0 ; i0< t->fiber_count[0]; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					TYPE* pr = partial_results + r;
					TYPE tval = t->val[i2];
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					//#pragma omp simd
					//printf("i vals are %d %d %d %d\n",i0,i1,i2,i3);
					//printf("%d %d %d %d\n",t->ind[0][i0],t->ind[1][i1],t->ind[2][i2],t->ind[3][i3]);
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						pr[y] += tval * matval[y];	// TTM step			
					}
				}	
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				TYPE* intval = t->intval[1] + i1*r;
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					//printf("1st level loop %lf\n",partial_results[r+y]);
					partial_results[y] += partial_results[r+y] * matval[y]; // TTV
					intval[y] = partial_results[r + y];
					//partial_results[r+y] = 0;
				}
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 

			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}

		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

	}
	rem(partial_results_all);	
	return 0;

}
int mttkrp_hardwired_first_4(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD;
	#ifdef OMP
	num_th = omp_get_max_threads();
	//printf("here \n");
	#endif

	printf("num ths %d\n", num_th);

	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

	for(int i = 0 ; i < num_th*partial_results_size ; i++)
		partial_results_all[i] = 0;

	set_matrix(*mats[0],0);

	#ifdef OMP
	#pragma omp parallel 
	#endif
	{
		int th = 0;
		#ifdef OMP
		th = omp_get_thread_num();
		#endif

		TYPE* partial_results = partial_results_all + th*partial_results_size;

		auto time_start = std::chrono::high_resolution_clock::now();
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0 ; i0< t->fiber_count[0]; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						TYPE* pr = partial_results + 2*r;
						TYPE tval = t->val[i3];
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						//#pragma omp simd
						//printf("i vals are %d %d %d %d\n",i0,i1,i2,i3);
						//printf("%d %d %d %d\n",t->ind[0][i0],t->ind[1][i1],t->ind[2][i2],t->ind[3][i3]);
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							pr[y] += tval * matval[y];	// TTM step			
						}
					}
					// write to intval
					
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					TYPE* intval = t->intval[2] + i2*r;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						//printf("2nd level loop %lf\n",partial_results[2*r+y]);
						partial_results[r+y] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y]= partial_results[2*r + y];
						
					}
					#pragma omp simd
					for(int y=0; y<r; y++)
					{
						partial_results[2*r+y] = 0;
					}
				}
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				TYPE* intval = t->intval[1] + i1*r;
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					//printf("1st level loop %lf\n",partial_results[r+y]);
					partial_results[y] += partial_results[r+y] * matval[y]; // TTV
					intval[y] = partial_results[r + y];
					//partial_results[r+y] = 0;
				}
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 

			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}

		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

	}
	rem(partial_results_all);	
	return 0;

}
int mttkrp_hardwired_first_5(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD;
	#ifdef OMP
	num_th = omp_get_max_threads();
	//printf("here \n");
	#endif

	printf("num ths %d\n", num_th);

	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

	for(int i = 0 ; i < num_th*partial_results_size ; i++)
		partial_results_all[i] = 0;

	set_matrix(*mats[0],0);

	#ifdef OMP
	#pragma omp parallel 
	#endif
	{
		int th = 0;
		#ifdef OMP
		th = omp_get_thread_num();
		#endif

		TYPE* partial_results = partial_results_all + th*partial_results_size;

		auto time_start = std::chrono::high_resolution_clock::now();
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0 ; i0< t->fiber_count[0]; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							TYPE* pr = partial_results + 3*r;
							TYPE tval = t->val[i4];
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							//#pragma omp simd
							//printf("i vals are %d %d %d %d\n",i0,i1,i2,i3);
							//printf("%d %d %d %d\n",t->ind[0][i0],t->ind[1][i1],t->ind[2][i2],t->ind[3][i3]);
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr[y] += tval * matval[y];	// TTM step			
							}
						}
						// write to intval
						
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							//printf("2nd level loop %lf\n",partial_results[2*r+y]);
							partial_results[2*r+y] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y]= partial_results[3*r + y];
							
						}
						#pragma omp simd
						for(int y=0; y<r; y++)
						{
							partial_results[3*r+y] = 0;
						}


					}
					// write to intval
					
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					TYPE* intval = t->intval[2] + i2*r;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						//printf("2nd level loop %lf\n",partial_results[2*r+y]);
						partial_results[r+y] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y]= partial_results[2*r + y];
						
					}
					#pragma omp simd
					for(int y=0; y<r; y++)
					{
						partial_results[2*r+y] = 0;
					}
				}
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				TYPE* intval = t->intval[1] + i1*r;
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					//printf("1st level loop %lf\n",partial_results[r+y]);
					partial_results[y] += partial_results[r+y] * matval[y]; // TTV
					intval[y] = partial_results[r + y];
					//partial_results[r+y] = 0;
				}
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 

			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}

		}

		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

	}
	rem(partial_results_all);	
	return 0;

}


int mttkrp_hardwired_last_3(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all;

	int nmode;
	int num_th;
	
	

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();
	#else
		num_th = 1;
		int th = 0;
	#endif

	int partial_products_size = nmode*r + PAD;
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

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
			if(VERBOSE == VERBOSE_HIGH)
				printf("th id is %d\n",th);
		#endif


		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		//TYPE* temp_res = temp_res_all + th*(r+64);
		
		TYPE* vals;
		if(t->num_th > 1)
		{
			vals = t->private_mats[th]->val;
		}
		else
		{
			vals = mats[mode]->val;
		}

		memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

	
		auto time_start = std::chrono::high_resolution_clock::now();
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0; i0 < (t->fiber_count[0]) ; i0++)
		{
			// init first pp
			//printf("first index is %d\n",i0);
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);

			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//printf(" middle index is %d\n",i1);

				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+r] = partial_products[y] * matval1[y];	
				}

				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1]  ; i2++)
				{
					const idx_t row_id = t->ind[2][i2];
					//printf("  last index is %d stopping at %d\n",i2,t->ptr[1][i1+1]);
					TYPE* xx  = vals + t->ind[2][i2]*(mats[mode]->dim2);
					TYPE* yy = partial_products + r ;

					TYPE tval = t->val[i2];


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
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}

	LIKWID_MARKER_CLOSE;

	t->num_th = num_th;
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}


	rem(partial_products_all);
	

	return 0;
}

int mttkrp_hardwired_last_4(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all;

	int nmode;
	int num_th;
	
	

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();
	#else
		num_th = 1;
		int th = 0;
	#endif

	int partial_products_size = nmode*r + PAD;
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

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
			if(VERBOSE == VERBOSE_HIGH)
				printf("th id is %d\n",th);
		#endif


		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		//TYPE* temp_res = temp_res_all + th*(r+64);
		
		TYPE* vals;
		if(t->num_th > 1)
		{
			vals = t->private_mats[th]->val;
		}
		else
		{
			vals = mats[mode]->val;
		}

		memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

		auto time_start = std::chrono::high_resolution_clock::now();
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0; i0 < (t->fiber_count[0]) ; i0++)
		{
			// init first pp
			//printf("first index is %d\n",i0);
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);

			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//printf(" middle index is %d\n",i1);

				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+r] = partial_products[y] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1]  ; i2++)
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//printf(" middle index is %d\n",i1);

					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+r] * matval2[y];	
					}

					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1]  ; i3++)
					{
						const idx_t row_id = t->ind[3][i3];
						//printf("  last index is %d stopping at %d\n",i2,t->ptr[1][i1+1]);
						TYPE* xx  = vals + row_id*(mats[mode]->dim2);
						TYPE* yy = partial_products + 2*r ;

						TYPE tval = t->val[i3];
						//#ifdef OMP
						//mutex_set_lock(mutex,row_id);
						//#endif

						#pragma omp simd
						for(int i=0 ; i<r ; i++)
						{
							TYPE increment = yy [i] * tval;
							xx [i]	+= increment;
							//#pragma omp critical
							//printf("writing to pos %llu and row %d in thread %d val is %lf\n",xx,row_id,th,increment);
							//printf("last mode %lf %lf %lf to pos mat[%d][%d][%d] \n", xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
						}
						//#ifdef OMP
						//mutex_unset_lock(mutex,row_id);
						//#endif
					}
				}
			}
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}

	
	LIKWID_MARKER_CLOSE;

	t->num_th = num_th;
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}


	rem(partial_products_all);
	

	return 0;
}

int mttkrp_hardwired_last_5(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all;

	int nmode;
	int num_th;
	
	

	nmode = t->nmode;
	#ifdef OMP
		num_th = omp_get_max_threads();
	#else
		num_th = 1;
		int th = 0;
	#endif

	int partial_products_size = nmode*r + PAD;
	/*	
	for(i=0 ;i<nmode; i++)
	{
		print_matrix(*mats[i]);
	}
	*/
	printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	//vals = (TYPE* ) malloc(mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	//vals = mats[mode]->val;
	//for(int i=0 ; i<(mats[mode]->dim1)*(mats[mode]->dim2) ; i++)
	//	vals[i] = 0;
	
	//memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

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
			if(VERBOSE == VERBOSE_HIGH)
				printf("th id is %d\n",th);
		#endif


		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		//TYPE* temp_res = temp_res_all + th*(r+64);
		
		TYPE* vals;
		if(t->num_th > 1)
		{
			vals = t->private_mats[th]->val;
		}
		else
		{
			vals = mats[mode]->val;
		}

		memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));


		auto time_start = std::chrono::high_resolution_clock::now();
				
		#ifdef OMP
		#pragma omp for schedule(dynamic,1)
		#endif
		for(idx_t i0 = 0; i0 < (t->fiber_count[0]) ; i0++)
		{
			// init first pp
			//printf("first index is %d\n",i0);
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);

			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//printf(" middle index is %d\n",i1);

				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+r] = partial_products[y] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1]  ; i2++)
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//printf(" middle index is %d\n",i1);

					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+r] * matval2[y];	
					}

					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1]  ; i3++)
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//printf(" middle index is %d\n",i1);

						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}

						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1]  ; i4++)
						{
							const idx_t row_id = t->ind[4][i4];
							//printf("  last index is %d stopping at %d\n",i2,t->ptr[1][i1+1]);
							TYPE* xx  = vals + row_id*(mats[mode]->dim2);
							TYPE* yy = partial_products + 3*r ;

							TYPE tval = t->val[i4];


							#pragma omp simd
							for(int i=0 ; i<r ; i++)
							{
								TYPE increment = yy [i] * tval;
								xx [i]	+= increment;
								//printf("last mode %lf %lf %lf to pos mat[%d][%d][%d] \n", xx [i], yy[i] , tval, mode, t->ind[nmode-1][it],i);
							}

						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}

	LIKWID_MARKER_CLOSE;

	t->num_th = num_th;
	if(num_th > 1)
	{
		auto time_start = std::chrono::high_resolution_clock::now();
		reduce(t,r,mats[mode]);
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		


	}

	rem(partial_products_all);
	

	return 0;
}

#endif