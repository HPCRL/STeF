#ifndef MTTKRP_COMBINED_CPP
#define MTTKRP_COMBINED_CPP

#include "../inc/mttkrp_combined.h"

int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile , int mode, bool intv1 )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; idx_t* thread_start = t->thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

	// Initialize the result arrays to 0
	if (mode > 0)
	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	else
		set_matrix(*mats[0],0);
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
		partial_results_all[i] = 0;
	
	
	
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
			TYPE* mv1 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* pr;
				TYPE* mv2 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if (!intv1 || mode == 2)
				{
					pr = partial_results;
					memset(pr,0,sizeof(TYPE)*r);
				}
				else
					pr = t->intval[1] + i1*r;

				if(mode == 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr[y] += mv1[y] * mv2[y];	// KrP step			
					}
				}

				if( mode != 1 || !intv1 )
				{
					
					for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
					{	
						TYPE tval = t->val[i2];
						
						
						if(mode != 2)
						{
							TYPE* matval = mats[2]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr[y] += tval * matval[y];	// dot TTM step			
							}
						}
						if(mode == 2)
						{
							TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								matval[y] += tval * pr[y];	// saxpy step			
							}
						}
					}	
				}			
				
				if(mode == 0)
				{
					TYPE* matval = mv1;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						mv1[y] += pr[y] * mv2[y]; // dot TTV
					}
				}
				if (mode == 1)
				{
					TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						matval[y] += pr[y] * mv1[y]; // saxpy
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
	
	if (mode > 0)
	{
		reduce(t,r,mats[mode]);
	}
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}

/*
template <int mode, bool intv1, bool intv2>
int mttkrp_combined_4(csf* t, int r, matrix** mats, int profile )
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


template <int mode, bool intv1, bool intv2, bool intv3>
int mttkrp_combined_5(csf* t, int r, matrix** mats, int profile )
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
		std:IKWID_MARKER_STOP("Compute");
		}
	}
	
	
	
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}
*/

int mttkrp_combined(csf* t, int r, matrix** mats, int profile , int mode, bool privatized, bool* intv )
{
	if(t->nmode == 3)
		if(mode == 0)
			if(intv[0])			
				if(privatized)
					mttkrp_combined_3<0,true,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<0,true,false>(t,r,mats,profile);
			else
				if(privatized)
					mttkrp_combined_3<0,false,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<0,false,false>(t,r,mats,profile);
		else if (mode == 1)
			if(intv[0])			
				if(privatized)
					mttkrp_combined_3<1,true,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<1,true,false>(t,r,mats,profile);
			else
				if(privatized)
					mttkrp_combined_3<1,false,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<1,false,false>(t,r,mats,profile);
		else 
			if(intv[0])			
				if(privatized)
					mttkrp_combined_3<2,true,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<2,true,false>(t,r,mats,profile);
			else
				if(privatized)
					mttkrp_combined_3<2,false,true>(t,r,mats,profile);
				else
					mttkrp_combined_3<2,false,false>(t,r,mats,profile);

	return 0;
}

#endif