#ifndef MTTKRP_HARDWIRED_CPP
#define MTTKRP_HARDWIRED_CPP
#include "../inc/mttkrp_hardwired.h"

#define SCHEDULE dynamic
#define CHUNK_SIZE 16




void find_thread_start(csf* t)
{
	idx_t* thread_start = t->thread_start;
	int nmode = t->nmode;
	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	idx_t total_work = 0;

	idx_t nnz = t->fiber_count[nmode-1];
	thread_start = (idx_t*) malloc(sizeof(idx_t)*(num_th+1));

	memset(thread_start,0,sizeof(idx_t)*(num_th+1));

	for ( int d = 0; d < nmode ; d++)
	{
		total_work += t->fiber_count[d];
	}

	idx_t sid = 0;
	idx_t partial_work = 0;
	idx_t max_slice_work = 0;

	for(int th = 0 ; th < num_th ; th++)
	{
		//idx_t partial_work = 0;
		while(sid < t->fiber_count[0] && partial_work < (total_work/ num_th) * (th+1) )
		{
			idx_t loc1 = sid, loc2 = sid + 1;
			idx_t local_work = 0;
			for ( int d = 0; d < nmode ; d++)
			{
				partial_work += loc2 - loc1;
				local_work += loc2 - loc1;

				if(d < nmode-1)
				{
//					partial_work += loc2 - loc1;
					loc1 = t->ptr[d][loc1];
					loc2 = t->ptr[d][loc2];
					
				}	
				if(local_work > max_slice_work)
					max_slice_work = local_work;
			}	
			sid ++;
		}

		thread_start[th+1] = sid;
	}
	t->thread_start = thread_start;
	thread_start[num_th] = t->fiber_count[0];

	printf("Max work per slice is %lld, total work is %lld, avg work per thread is %lld \n",max_slice_work,total_work,total_work/num_th);

	for(int i=0 ; i<= num_th ; i++)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("ths %d %lld\n  ",i,thread_start[i]);
	}
}


int reduce(csf* t, int r, matrix* mat)
{
	// Reset the values in the result matrix to 0
	// memset (mat->val, 0 , sizeof(TYPE)*(mat->dim1)*(mat->dim2));
	int parallel_threshold = 1000;
	if (mat->dim1 > parallel_threshold)
	{
		#ifdef OMP
		#pragma omp parallel for
		#endif
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
					//if(VERBOSE >= VERBOSE_DEBUG) printf("reducing %lf %lf in th %d\n", outval[y], reduceval[y],i);
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
					//if(VERBOSE >= VERBOSE_DEBUG) printf("reducing %lf %lf in th %d\n", outval[y], reduceval[y],i);
					//reduceval[y] = 0;
				}
			}	
		}
	}


	return 0;
}
int mttkrp_hardwired_first_2(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* pr = partial_results + 0 * r;
				TYPE tval = t->val[i1];
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					pr[y] += tval * matval[y];	// TTM step			
				}
			}	
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_3(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			TYPE* mv = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* pr = t->intval[1] + i1*r;

				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{	
					TYPE tval = t->val[i2];
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					
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

					mv[y] += pr[y] * matval[y]; // TTV
				//	intval[y] = partial_results[1*r + y];
				}
				/*
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
				*/
			}
			// write to output matrix
			/*
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			*/
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_4(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start; 
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						TYPE* pr = partial_results + 2 * r;
						TYPE tval = t->val[i3];
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_5(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							TYPE* pr = partial_results + 3 * r;
							TYPE tval = t->val[i4];
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							
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
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_6(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								TYPE* pr = partial_results + 4 * r;
								TYPE tval = t->val[i5];
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									pr[y] += tval * matval[y];	// TTM step			
								}
							}	
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[y + 3*r] += partial_results[4*r+y] * matval[y]; // TTV
								intval[y] = partial_results[4*r + y];
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_7(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									TYPE* pr = partial_results + 5 * r;
									TYPE tval = t->val[i6];
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr[y] += tval * matval[y];	// TTM step			
									}
								}	
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[y + 4*r] += partial_results[5*r+y] * matval[y]; // TTV
									intval[y] = partial_results[5*r + y];
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[y + 3*r] += partial_results[4*r+y] * matval[y]; // TTV
								intval[y] = partial_results[4*r + y];
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_8(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										TYPE* pr = partial_results + 6 * r;
										TYPE tval = t->val[i7];
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											pr[y] += tval * matval[y];	// TTM step			
										}
									}	
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[y + 5*r] += partial_results[6*r+y] * matval[y]; // TTV
										intval[y] = partial_results[6*r + y];
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[y + 4*r] += partial_results[5*r+y] * matval[y]; // TTV
									intval[y] = partial_results[5*r + y];
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[y + 3*r] += partial_results[4*r+y] * matval[y]; // TTV
								intval[y] = partial_results[4*r + y];
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_9(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										for(idx_t i8 = t->ptr[7][i7] ; i8< t->ptr[7][i7+1]; i8++)
										{
											TYPE* pr = partial_results + 7 * r;
											TYPE tval = t->val[i8];
											TYPE* matval = (mats[8]->val) + ((mats[8]) -> dim2) * t->ind[8][i8];
											
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												pr[y] += tval * matval[y];	// TTM step			
											}
										}	
										// write to intval
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										TYPE* intval = t->intval[7] + i7*r;
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											partial_results[y + 6*r] += partial_results[7*r+y] * matval[y]; // TTV
											intval[y] = partial_results[7*r + y];
										}
										
										#pragma omp simd
										for(int y=0; y<r; y++)
										{
											partial_results[7*r+y] = 0;
										}
									}
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[y + 5*r] += partial_results[6*r+y] * matval[y]; // TTV
										intval[y] = partial_results[6*r + y];
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[y + 4*r] += partial_results[5*r+y] * matval[y]; // TTV
									intval[y] = partial_results[5*r + y];
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[y + 3*r] += partial_results[4*r+y] * matval[y]; // TTV
								intval[y] = partial_results[4*r + y];
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_10(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
	}
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										for(idx_t i8 = t->ptr[7][i7] ; i8< t->ptr[7][i7+1]; i8++)
										{
											for(idx_t i9 = t->ptr[8][i8] ; i9< t->ptr[8][i8+1]; i9++)
											{
												TYPE* pr = partial_results + 8 * r;
												TYPE tval = t->val[i9];
												TYPE* matval = (mats[9]->val) + ((mats[9]) -> dim2) * t->ind[9][i9];
												
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													pr[y] += tval * matval[y];	// TTM step			
												}
											}	
											// write to intval
											TYPE* matval = (mats[8]->val) + ((mats[8]) -> dim2) * t->ind[8][i8];
											TYPE* intval = t->intval[8] + i8*r;
											
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												partial_results[y + 7*r] += partial_results[8*r+y] * matval[y]; // TTV
												intval[y] = partial_results[8*r + y];
											}
											
											#pragma omp simd
											for(int y=0; y<r; y++)
											{
												partial_results[8*r+y] = 0;
											}
										}
										// write to intval
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										TYPE* intval = t->intval[7] + i7*r;
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											partial_results[y + 6*r] += partial_results[7*r+y] * matval[y]; // TTV
											intval[y] = partial_results[7*r + y];
										}
										
										#pragma omp simd
										for(int y=0; y<r; y++)
										{
											partial_results[7*r+y] = 0;
										}
									}
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[y + 5*r] += partial_results[6*r+y] * matval[y]; // TTV
										intval[y] = partial_results[6*r + y];
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[y + 4*r] += partial_results[5*r+y] * matval[y]; // TTV
									intval[y] = partial_results[5*r + y];
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[y + 3*r] += partial_results[4*r+y] * matval[y]; // TTV
								intval[y] = partial_results[4*r + y];
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[y + 2*r] += partial_results[3*r+y] * matval[y]; // TTV
							intval[y] = partial_results[3*r + y];
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
						partial_results[y + 1*r] += partial_results[2*r+y] * matval[y]; // TTV
						intval[y] = partial_results[2*r + y];
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
					partial_results[y + 0*r] += partial_results[1*r+y] * matval[y]; // TTV
					intval[y] = partial_results[1*r + y];
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_2(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* pr = partial_results + 0 * r;
				TYPE tval = t->val[i1];
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					pr[y] += tval * matval[y];	// TTM step			
				}
			}	
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_3(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					TYPE* pr = partial_results + 1 * r;
					TYPE tval = t->val[i2];
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						pr[y] += tval * matval[y];	// TTM step			
					}
				}	
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_4(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	partial_results_all[i] = 0;
	
	set_matrix(*mats[0],0);

	if(profile == mode)
        {
                if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		int th = 0;
		#ifdef OMP
		th = omp_get_thread_num();
		#endif
		auto time_start = std::chrono::high_resolution_clock::now();
		TYPE* partial_results = partial_results_all + th*partial_results_size;
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						TYPE* pr = partial_results + 2 * r;
						TYPE tval = t->val[i3];
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							pr[y] += tval * matval[y];	// TTM step			
						}
					}	
					// write to intval
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					
					
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
					}
					
					#pragma omp simd
					for(int y=0; y<r; y++)
					{
						partial_results[2*r+y] = 0;
					}
				}
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

		if(profile == mode)
                {
                        LIKWID_MARKER_STOP("Compute");
                }
        }

        LIKWID_MARKER_CLOSE; 

		
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_5(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							TYPE* pr = partial_results + 3 * r;
							TYPE tval = t->val[i4];
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr[y] += tval * matval[y];	// TTM step			
							}
						}	
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
						}
						
						#pragma omp simd
						for(int y=0; y<r; y++)
						{
							partial_results[3*r+y] = 0;
						}
					}
					// write to intval
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					
					
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
					}
					
					#pragma omp simd
					for(int y=0; y<r; y++)
					{
						partial_results[2*r+y] = 0;
					}
				}
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_6(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								TYPE* pr = partial_results + 4 * r;
								TYPE tval = t->val[i5];
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									pr[y] += tval * matval[y];	// TTM step			
								}
							}	
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[3*r + y] += partial_results[4*r+y] * matval[y]; // TTV
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
						}
						
						#pragma omp simd
						for(int y=0; y<r; y++)
						{
							partial_results[3*r+y] = 0;
						}
					}
					// write to intval
					TYPE* matval = (mats[2]->val) + ((mats[2]) -> dim2) * t->ind[2][i2];
					
					
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
					}
					
					#pragma omp simd
					for(int y=0; y<r; y++)
					{
						partial_results[2*r+y] = 0;
					}
				}
				// write to intval
				TYPE* matval = (mats[1]->val) + ((mats[1]) -> dim2) * t->ind[1][i1];
				
				
				#pragma omp simd
				for(int y=0 ; y<r ; y++)
				{
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_7(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									TYPE* pr = partial_results + 5 * r;
									TYPE tval = t->val[i6];
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr[y] += tval * matval[y];	// TTM step			
									}
								}	
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[4*r + y] += partial_results[5*r+y] * matval[y]; // TTV
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[3*r + y] += partial_results[4*r+y] * matval[y]; // TTV
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
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
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
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
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_8(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										TYPE* pr = partial_results + 6 * r;
										TYPE tval = t->val[i7];
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											pr[y] += tval * matval[y];	// TTM step			
										}
									}	
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[5*r + y] += partial_results[6*r+y] * matval[y]; // TTV
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[4*r + y] += partial_results[5*r+y] * matval[y]; // TTV
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[3*r + y] += partial_results[4*r+y] * matval[y]; // TTV
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
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
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
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
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_9(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										for(idx_t i8 = t->ptr[7][i7] ; i8< t->ptr[7][i7+1]; i8++)
										{
											TYPE* pr = partial_results + 7 * r;
											TYPE tval = t->val[i8];
											TYPE* matval = (mats[8]->val) + ((mats[8]) -> dim2) * t->ind[8][i8];
											
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												pr[y] += tval * matval[y];	// TTM step			
											}
										}	
										// write to intval
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										TYPE* intval = t->intval[7] + i7*r;
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											partial_results[6*r + y] += partial_results[7*r+y] * matval[y]; // TTV
										}
										
										#pragma omp simd
										for(int y=0; y<r; y++)
										{
											partial_results[7*r+y] = 0;
										}
									}
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[5*r + y] += partial_results[6*r+y] * matval[y]; // TTV
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[4*r + y] += partial_results[5*r+y] * matval[y]; // TTV
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[3*r + y] += partial_results[4*r+y] * matval[y]; // TTV
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
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
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
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
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_first_not_fused_10(csf* t, int mode, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
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
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
				{
					for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
					{
						for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
						{
							for(idx_t i5 = t->ptr[4][i4] ; i5< t->ptr[4][i4+1]; i5++)
							{
								for(idx_t i6 = t->ptr[5][i5] ; i6< t->ptr[5][i5+1]; i6++)
								{
									for(idx_t i7 = t->ptr[6][i6] ; i7< t->ptr[6][i6+1]; i7++)
									{
										for(idx_t i8 = t->ptr[7][i7] ; i8< t->ptr[7][i7+1]; i8++)
										{
											for(idx_t i9 = t->ptr[8][i8] ; i9< t->ptr[8][i8+1]; i9++)
											{
												TYPE* pr = partial_results + 8 * r;
												TYPE tval = t->val[i9];
												TYPE* matval = (mats[9]->val) + ((mats[9]) -> dim2) * t->ind[9][i9];
												
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													pr[y] += tval * matval[y];	// TTM step			
												}
											}	
											// write to intval
											TYPE* matval = (mats[8]->val) + ((mats[8]) -> dim2) * t->ind[8][i8];
											TYPE* intval = t->intval[8] + i8*r;
											
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												partial_results[7*r + y] += partial_results[8*r+y] * matval[y]; // TTV
											}
											
											#pragma omp simd
											for(int y=0; y<r; y++)
											{
												partial_results[8*r+y] = 0;
											}
										}
										// write to intval
										TYPE* matval = (mats[7]->val) + ((mats[7]) -> dim2) * t->ind[7][i7];
										TYPE* intval = t->intval[7] + i7*r;
										
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											partial_results[6*r + y] += partial_results[7*r+y] * matval[y]; // TTV
										}
										
										#pragma omp simd
										for(int y=0; y<r; y++)
										{
											partial_results[7*r+y] = 0;
										}
									}
									// write to intval
									TYPE* matval = (mats[6]->val) + ((mats[6]) -> dim2) * t->ind[6][i6];
									TYPE* intval = t->intval[6] + i6*r;
									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										partial_results[5*r + y] += partial_results[6*r+y] * matval[y]; // TTV
									}
									
									#pragma omp simd
									for(int y=0; y<r; y++)
									{
										partial_results[6*r+y] = 0;
									}
								}
								// write to intval
								TYPE* matval = (mats[5]->val) + ((mats[5]) -> dim2) * t->ind[5][i5];
								TYPE* intval = t->intval[5] + i5*r;
								
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									partial_results[4*r + y] += partial_results[5*r+y] * matval[y]; // TTV
								}
								
								#pragma omp simd
								for(int y=0; y<r; y++)
								{
									partial_results[5*r+y] = 0;
								}
							}
							// write to intval
							TYPE* matval = (mats[4]->val) + ((mats[4]) -> dim2) * t->ind[4][i4];
							TYPE* intval = t->intval[4] + i4*r;
							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								partial_results[3*r + y] += partial_results[4*r+y] * matval[y]; // TTV
							}
							
							#pragma omp simd
							for(int y=0; y<r; y++)
							{
								partial_results[4*r+y] = 0;
							}
						}
						// write to intval
						TYPE* matval = (mats[3]->val) + ((mats[3]) -> dim2) * t->ind[3][i3];
						TYPE* intval = t->intval[3] + i3*r;
						
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							partial_results[2*r + y] += partial_results[3*r+y] * matval[y]; // TTV
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
						partial_results[1*r + y] += partial_results[2*r+y] * matval[y]; // TTV
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
					partial_results[0*r + y] += partial_results[1*r+y] * matval[y]; // TTV
				}
				
				#pragma omp simd
				for(int y=0; y<r; y++)
				{
					partial_results[1*r+y] = 0;
				}
			}
			// write to output matrix
			TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 
			
			#pragma omp simd
			for(int y=0 ; y<r ; y++)
			{
				//	if(VERBOSE >= VERBOSE_DEBUG) printf("0th level loop %lf\n",partial_results[y]);
				matval[y] = partial_results[y];
				partial_results[y] = 0;
			}
			
		}
		
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
	}
	rem(partial_results_all);	
	return 0;
	
}
int mttkrp_hardwired_last_2(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1]  ; i1++)
			{
				const idx_t row_id = t->ind[1][i1];
				TYPE* xx  = vals + t->ind[1][i1]*(mats[mode]->dim2);
				TYPE* yy = partial_products + 0*r ;
				TYPE tval = t->val[i1];
				
				
				#pragma omp simd
				for(int i=0 ; i<r ; i++)
				{
					TYPE increment = yy [i] * tval;
					xx [i]	+= increment;
				}
				
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_3(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
//	TYPE const * const  partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	idx_t* thread_start = t->thread_start;
	int nmode;
	int num_th;
	nmode = t->nmode;

	//nst int r = 32;
	//#define r 32
	
	#ifdef OMP
	num_th = omp_get_max_threads();
	#else
	num_th = 1;
	const int th = 0;
	#endif
	
	int partial_products_size = nmode*r + PAD;
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	TYPE  * const __restrict__  partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
		LIKWID_MARKER_INIT;
	}

//	memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	
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
//		th = omp_get_thread_num();
//                if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
	}

	TYPE const * const __restrict__ mat0 = mats[0]->val;	
	TYPE const * const __restrict__ mat1 = mats[1]->val;
	const idx_t dim2 =  mats[0]->dim2;
	idx_t const * const __restrict__ ind0 = (t->ind)[0];
        idx_t const * const __restrict__ ind1 = (t->ind)[1];
	idx_t const * const __restrict__ ind2 = (t->ind)[2];
	idx_t const * const __restrict__ ptr0 = (t->ptr)[0];
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
		#ifdef OMP
		const int th = omp_get_thread_num();
		if(VERBOSE == VERBOSE_HIGH)
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
//		TYPE * const __restrict__ partial_products = partial_products_all + th*partial_products_size;
//		TYPE * const __restrict__ partial_products = (TYPE* ) malloc((partial_products_size)*sizeof(TYPE));		
		TYPE * const __restrict__ pp_out = (TYPE* ) malloc((r)*sizeof(TYPE));
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
		const idx_t i0_start = thread_start[th];
		const idx_t i0_stop = thread_start[th+1];
		for(idx_t i0 = i0_start ; i0 < i0_stop ; i0++)
		{
			
			TYPE const * const __restrict__  matval0 = mat0 + (dim2 * ind0[i0]);
			
//			#pragma omp simd
//			for(int y=0; y<r ; y++)
//			{
//				partial_products[y] = matval0[y];	
//			}
			
			const idx_t i1_start = ptr0[i0];
			const idx_t i1_stop = ptr0[i0+1];			

			for(idx_t i1 = i1_start ; i1 < i1_stop ; i1++)	
			{
				TYPE const * const __restrict__  matval1 = mat1 + (dim2 * ind1[i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);

//				TYPE * const __restrict__ pp_out = partial_products + r;
				//TYPE const * const __restrict__ in = partial_products;

				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
//					pp_out[y] = in[y] * matval1[y];
					pp_out[y] = matval0[y] * matval1[y];					
				}
				
//				#pragma omp simd
//				for(int y=0; y<r ; y++)
//				{
//					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
//				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1]  ; i2++)
				{
					const idx_t row_id = t->ind[2][i2];
					TYPE* const __restrict__ xx  = vals + ind2[i2]*dim2;
					//TYPE const * const __restrict__ yy = partial_products + 1*r ;
					TYPE tval = t->val[i2];
					
					
					#pragma omp simd
					for(int i=0 ; i<r ; i++)
					{
						TYPE increment = pp_out [i] * tval;
						xx [i]	+= increment;
					}
					
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
//		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
//	rem(partial_products_all);
	free(partial_products_all);	
	return 0;
}
int mttkrp_hardwired_last_4(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1]  ; i3++)
					{
						const idx_t row_id = t->ind[3][i3];
						TYPE* xx  = vals + t->ind[3][i3]*(mats[mode]->dim2);
						TYPE* yy = partial_products + 2*r ;
						TYPE tval = t->val[i3];
						
						
						#pragma omp simd
						for(int i=0 ; i<r ; i++)
						{
							TYPE increment = yy [i] * tval;
							xx [i]	+= increment;
						}
						
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_5(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1]  ; i4++)
						{
							const idx_t row_id = t->ind[4][i4];
							TYPE* xx  = vals + t->ind[4][i4]*(mats[mode]->dim2);
							TYPE* yy = partial_products + 3*r ;
							TYPE tval = t->val[i4];
							
							
							#pragma omp simd
							for(int i=0 ; i<r ; i++)
							{
								TYPE increment = yy [i] * tval;
								xx [i]	+= increment;
							}
							
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_6(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1]  ; i5++)
							{
								const idx_t row_id = t->ind[5][i5];
								TYPE* xx  = vals + t->ind[5][i5]*(mats[mode]->dim2);
								TYPE* yy = partial_products + 4*r ;
								TYPE tval = t->val[i5];
								
								
								#pragma omp simd
								for(int i=0 ; i<r ; i++)
								{
									TYPE increment = yy [i] * tval;
									xx [i]	+= increment;
								}
								
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_7(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1]  ; i6++)
								{
									const idx_t row_id = t->ind[6][i6];
									TYPE* xx  = vals + t->ind[6][i6]*(mats[mode]->dim2);
									TYPE* yy = partial_products + 5*r ;
									TYPE tval = t->val[i6];
									
									
									#pragma omp simd
									for(int i=0 ; i<r ; i++)
									{
										TYPE increment = yy [i] * tval;
										xx [i]	+= increment;
									}
									
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_8(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1]  ; i7++)
									{
										const idx_t row_id = t->ind[7][i7];
										TYPE* xx  = vals + t->ind[7][i7]*(mats[mode]->dim2);
										TYPE* yy = partial_products + 6*r ;
										TYPE tval = t->val[i7];
										
										
										#pragma omp simd
										for(int i=0 ; i<r ; i++)
										{
											TYPE increment = yy [i] * tval;
											xx [i]	+= increment;
										}
										
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_9(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1] ; i7++)	
									{
										TYPE* matval7 = mats[7]->val + ((mats[7]->dim2) * t->ind[7][i7]);
										//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
										
										#pragma omp simd
										for(int y=0; y<r ; y++)
										{
											partial_products[y+7*r] = partial_products[y+6*r] * matval7[y];	
										}
										for(idx_t i8 = t->ptr[7][i7]; i8 < t->ptr[7][i7+1]  ; i8++)
										{
											const idx_t row_id = t->ind[8][i8];
											TYPE* xx  = vals + t->ind[8][i8]*(mats[mode]->dim2);
											TYPE* yy = partial_products + 7*r ;
											TYPE tval = t->val[i8];
											
											
											#pragma omp simd
											for(int i=0 ; i<r ; i++)
											{
												TYPE increment = yy [i] * tval;
												xx [i]	+= increment;
											}
											
										}
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_10(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1] ; i7++)	
									{
										TYPE* matval7 = mats[7]->val + ((mats[7]->dim2) * t->ind[7][i7]);
										//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
										
										#pragma omp simd
										for(int y=0; y<r ; y++)
										{
											partial_products[y+7*r] = partial_products[y+6*r] * matval7[y];	
										}
										for(idx_t i8 = t->ptr[7][i7]; i8 < t->ptr[7][i7+1] ; i8++)	
										{
											TYPE* matval8 = mats[8]->val + ((mats[8]->dim2) * t->ind[8][i8]);
											//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
											
											#pragma omp simd
											for(int y=0; y<r ; y++)
											{
												partial_products[y+8*r] = partial_products[y+7*r] * matval8[y];	
											}
											for(idx_t i9 = t->ptr[8][i8]; i9 < t->ptr[8][i8+1]  ; i9++)
											{
												const idx_t row_id = t->ind[9][i9];
												TYPE* xx  = vals + t->ind[9][i9]*(mats[mode]->dim2);
												TYPE* yy = partial_products + 8*r ;
												TYPE tval = t->val[i9];
												
												
												#pragma omp simd
												for(int i=0 ; i<r ; i++)
												{
													TYPE increment = yy [i] * tval;
													xx [i]	+= increment;
												}
												
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_2(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1]  ; i1++)
			{
				const idx_t row_id = t->ind[1][i1];
				TYPE* xx  = vals + t->ind[1][i1]*(mats[mode]->dim2);
				TYPE* yy = partial_products + 0*r ;
				TYPE* intval = t->intval[1] + i1*r ;
				
				
				#pragma omp simd
				for(int i=0 ; i<r ; i++)
				{
					TYPE increment = yy [i] * intval[i];
					xx [i]	+= increment;
				}
				
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_3(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1]  ; i2++)
				{
					const idx_t row_id = t->ind[2][i2];
					TYPE* xx  = vals + t->ind[2][i2]*(mats[mode]->dim2);
					TYPE* yy = partial_products + 1*r ;
					TYPE* intval = t->intval[2] + i2*r ;
					
					
					#pragma omp simd
					for(int i=0 ; i<r ; i++)
					{
						TYPE increment = yy [i] * intval[i];
						xx [i]	+= increment;
					}
					
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_4(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1]  ; i3++)
					{
						const idx_t row_id = t->ind[3][i3];
						TYPE* xx  = vals + t->ind[3][i3]*(mats[mode]->dim2);
						TYPE* yy = partial_products + 2*r ;
						TYPE* intval = t->intval[3] + i3*r ;
						
						
						#pragma omp simd
						for(int i=0 ; i<r ; i++)
						{
							TYPE increment = yy [i] * intval[i];
							xx [i]	+= increment;
						}
						
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_5(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1]  ; i4++)
						{
							const idx_t row_id = t->ind[4][i4];
							TYPE* xx  = vals + t->ind[4][i4]*(mats[mode]->dim2);
							TYPE* yy = partial_products + 3*r ;
							TYPE* intval = t->intval[4] + i4*r ;
							
							
							#pragma omp simd
							for(int i=0 ; i<r ; i++)
							{
								TYPE increment = yy [i] * intval[i];
								xx [i]	+= increment;
							}
							
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_6(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1]  ; i5++)
							{
								const idx_t row_id = t->ind[5][i5];
								TYPE* xx  = vals + t->ind[5][i5]*(mats[mode]->dim2);
								TYPE* yy = partial_products + 4*r ;
								TYPE* intval = t->intval[5] + i5*r ;
								
								
								#pragma omp simd
								for(int i=0 ; i<r ; i++)
								{
									TYPE increment = yy [i] * intval[i];
									xx [i]	+= increment;
								}
								
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_7(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1]  ; i6++)
								{
									const idx_t row_id = t->ind[6][i6];
									TYPE* xx  = vals + t->ind[6][i6]*(mats[mode]->dim2);
									TYPE* yy = partial_products + 5*r ;
									TYPE* intval = t->intval[6] + i6*r ;
									
									
									#pragma omp simd
									for(int i=0 ; i<r ; i++)
									{
										TYPE increment = yy [i] * intval[i];
										xx [i]	+= increment;
									}
									
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_8(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1]  ; i7++)
									{
										const idx_t row_id = t->ind[7][i7];
										TYPE* xx  = vals + t->ind[7][i7]*(mats[mode]->dim2);
										TYPE* yy = partial_products + 6*r ;
										TYPE* intval = t->intval[7] + i7*r ;
										
										
										#pragma omp simd
										for(int i=0 ; i<r ; i++)
										{
											TYPE increment = yy [i] * intval[i];
											xx [i]	+= increment;
										}
										
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_9(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1] ; i7++)	
									{
										TYPE* matval7 = mats[7]->val + ((mats[7]->dim2) * t->ind[7][i7]);
										//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
										
										#pragma omp simd
										for(int y=0; y<r ; y++)
										{
											partial_products[y+7*r] = partial_products[y+6*r] * matval7[y];	
										}
										for(idx_t i8 = t->ptr[7][i7]; i8 < t->ptr[7][i7+1]  ; i8++)
										{
											const idx_t row_id = t->ind[8][i8];
											TYPE* xx  = vals + t->ind[8][i8]*(mats[mode]->dim2);
											TYPE* yy = partial_products + 7*r ;
											TYPE* intval = t->intval[8] + i8*r ;
											
											
											#pragma omp simd
											for(int i=0 ; i<r ; i++)
											{
												TYPE increment = yy [i] * intval[i];
												xx [i]	+= increment;
											}
											
										}
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
int mttkrp_hardwired_last_vec_10(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
	TYPE* partial_products_all; idx_t* thread_start = t->thread_start;
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
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));
	
	if(profile == mode)
	{
		if(VERBOSE >= VERBOSE_DEBUG) printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("th id is %d\n",th);
		#endif
		
		TYPE* partial_products;	
		partial_products = partial_products_all + th*partial_products_size;
		
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
				for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			
			TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);
			
			#pragma omp simd
			for(int y=0; y<r ; y++)
			{
				partial_products[y] = matval0[y];	
			}
			for(idx_t i1 = t->ptr[0][i0]; i1 < t->ptr[0][i0+1] ; i1++)	
			{
				TYPE* matval1 = mats[1]->val + ((mats[1]->dim2) * t->ind[1][i1]);
				//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
				
				#pragma omp simd
				for(int y=0; y<r ; y++)
				{
					partial_products[y+1*r] = partial_products[y+0*r] * matval1[y];	
				}
				for(idx_t i2 = t->ptr[1][i1]; i2 < t->ptr[1][i1+1] ; i2++)	
				{
					TYPE* matval2 = mats[2]->val + ((mats[2]->dim2) * t->ind[2][i2]);
					//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
					
					#pragma omp simd
					for(int y=0; y<r ; y++)
					{
						partial_products[y+2*r] = partial_products[y+1*r] * matval2[y];	
					}
					for(idx_t i3 = t->ptr[2][i2]; i3 < t->ptr[2][i2+1] ; i3++)	
					{
						TYPE* matval3 = mats[3]->val + ((mats[3]->dim2) * t->ind[3][i3]);
						//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
						
						#pragma omp simd
						for(int y=0; y<r ; y++)
						{
							partial_products[y+3*r] = partial_products[y+2*r] * matval3[y];	
						}
						for(idx_t i4 = t->ptr[3][i3]; i4 < t->ptr[3][i3+1] ; i4++)	
						{
							TYPE* matval4 = mats[4]->val + ((mats[4]->dim2) * t->ind[4][i4]);
							//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
							
							#pragma omp simd
							for(int y=0; y<r ; y++)
							{
								partial_products[y+4*r] = partial_products[y+3*r] * matval4[y];	
							}
							for(idx_t i5 = t->ptr[4][i4]; i5 < t->ptr[4][i4+1] ; i5++)	
							{
								TYPE* matval5 = mats[5]->val + ((mats[5]->dim2) * t->ind[5][i5]);
								//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
								
								#pragma omp simd
								for(int y=0; y<r ; y++)
								{
									partial_products[y+5*r] = partial_products[y+4*r] * matval5[y];	
								}
								for(idx_t i6 = t->ptr[5][i5]; i6 < t->ptr[5][i5+1] ; i6++)	
								{
									TYPE* matval6 = mats[6]->val + ((mats[6]->dim2) * t->ind[6][i6]);
									//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
									
									#pragma omp simd
									for(int y=0; y<r ; y++)
									{
										partial_products[y+6*r] = partial_products[y+5*r] * matval6[y];	
									}
									for(idx_t i7 = t->ptr[6][i6]; i7 < t->ptr[6][i6+1] ; i7++)	
									{
										TYPE* matval7 = mats[7]->val + ((mats[7]->dim2) * t->ind[7][i7]);
										//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
										
										#pragma omp simd
										for(int y=0; y<r ; y++)
										{
											partial_products[y+7*r] = partial_products[y+6*r] * matval7[y];	
										}
										for(idx_t i8 = t->ptr[7][i7]; i8 < t->ptr[7][i7+1] ; i8++)	
										{
											TYPE* matval8 = mats[8]->val + ((mats[8]->dim2) * t->ind[8][i8]);
											//if(VERBOSE >= VERBOSE_DEBUG) printf(" middle index is %d\n",i1);
											
											#pragma omp simd
											for(int y=0; y<r ; y++)
											{
												partial_products[y+8*r] = partial_products[y+7*r] * matval8[y];	
											}
											for(idx_t i9 = t->ptr[8][i8]; i9 < t->ptr[8][i8+1]  ; i9++)
											{
												const idx_t row_id = t->ind[9][i9];
												TYPE* xx  = vals + t->ind[9][i9]*(mats[mode]->dim2);
												TYPE* yy = partial_products + 8*r ;
												TYPE* intval = t->intval[9] + i9*r ;
												
												
												#pragma omp simd
												for(int i=0 ; i<r ; i++)
												{
													TYPE increment = yy [i] * intval[i];
													xx [i]	+= increment;
												}
												
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		
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
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
	}
	
	rem(partial_products_all);
	return 0;
}
#endif
