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
	{		
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
	}
	if(t->nmode == 4)
	{
		if( mode == 0 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<0, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<0, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<0, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_4<0, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<0, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<0, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<0, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && privatized)
			mttkrp_combined_4<0, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<1, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<1, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<1, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_4<1, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<1, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<1, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<1, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && privatized)
			mttkrp_combined_4<1, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<2, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<2, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<2, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_4<2, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<2, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<2, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<2, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && privatized)
			mttkrp_combined_4<2, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<3, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<3, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<3, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_4<3, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_4<3, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_4<3, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_4<3, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && privatized)
			mttkrp_combined_4<3, true, true, true>(t,r,mats,profile);
	}
	if(t->nmode == 5)
	{
		if( mode == 0 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<0, false, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<0, false, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<0, false, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<0, false, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<0, false, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<0, false, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<0, false, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<0, false, true, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<0, true, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<0, true, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<0, true, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<0, true, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<0, true, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<0, true, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<0, true, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<0, true, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<1, false, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<1, false, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<1, false, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<1, false, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<1, false, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<1, false, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<1, false, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<1, false, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<1, true, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<1, true, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<1, true, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<1, true, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<1, true, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<1, true, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<1, true, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<1, true, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<2, false, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<2, false, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<2, false, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<2, false, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<2, false, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<2, false, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<2, false, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<2, false, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<2, true, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<2, true, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<2, true, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<2, true, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<2, true, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<2, true, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<2, true, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<2, true, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<3, false, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<3, false, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<3, false, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<3, false, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<3, false, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<3, false, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<3, false, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<3, false, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<3, true, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<3, true, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<3, true, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<3, true, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<3, true, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<3, true, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<3, true, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<3, true, true, true, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<4, false, false, false, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<4, false, false, false, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<4, false, false, true, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<4, false, false, true, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<4, false, true, false, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<4, false, true, false, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<4, false, true, true, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<4, false, true, true, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<4, true, false, false, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<4, true, false, false, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<4, true, false, true, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_5<4, true, false, true, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_5<4, true, true, false, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_5<4, true, true, false, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_5<4, true, true, true, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_5<4, true, true, true, true>(t,r,mats,profile);
	}

	return 0;
}

int mttkrp_combined_lb(csf* t, int r, matrix** mats, int profile , int mode, bool privatized, bool* intv )
{
	if(t->nmode == 3)
	{
		if( mode == 0 && !intv[0] && !privatized)
			mttkrp_combined_lb_3<0, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && privatized)
			mttkrp_combined_lb_3<0, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !privatized)
			mttkrp_combined_lb_3<0, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && privatized)
			mttkrp_combined_lb_3<0, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !privatized)
			mttkrp_combined_lb_3<1, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && privatized)
			mttkrp_combined_lb_3<1, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !privatized)
			mttkrp_combined_lb_3<1, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && privatized)
			mttkrp_combined_lb_3<1, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !privatized)
			mttkrp_combined_lb_3<2, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && privatized)
			mttkrp_combined_lb_3<2, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !privatized)
			mttkrp_combined_lb_3<2, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && privatized)
			mttkrp_combined_lb_3<2, true, true>(t,r,mats,profile);
	}
	if(t->nmode == 4)
	{
		if( mode == 0 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<0, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<0, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<0, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<0, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<0, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<0, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<0, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<0, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<1, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<1, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<1, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<1, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<1, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<1, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<1, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<1, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<2, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<2, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<2, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<2, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<2, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<2, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<2, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<2, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<3, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<3, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<3, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<3, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !privatized)
			mttkrp_combined_lb_4<3, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && privatized)
			mttkrp_combined_lb_4<3, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !privatized)
			mttkrp_combined_lb_4<3, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && privatized)
			mttkrp_combined_lb_4<3, true, true, true>(t,r,mats,profile);
	}
	if(t->nmode == 5)
	{
		if( mode == 0 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<0, false, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<0, false, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<0, false, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<0, false, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<0, false, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<0, false, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<0, false, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<0, false, true, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<0, true, false, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<0, true, false, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<0, true, false, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<0, true, false, true, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<0, true, true, false, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<0, true, true, false, true>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<0, true, true, true, false>(t,r,mats,profile);
		else if( mode == 0 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<0, true, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<1, false, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<1, false, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<1, false, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<1, false, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<1, false, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<1, false, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<1, false, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<1, false, true, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<1, true, false, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<1, true, false, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<1, true, false, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<1, true, false, true, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<1, true, true, false, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<1, true, true, false, true>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<1, true, true, true, false>(t,r,mats,profile);
		else if( mode == 1 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<1, true, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<2, false, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<2, false, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<2, false, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<2, false, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<2, false, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<2, false, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<2, false, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<2, false, true, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<2, true, false, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<2, true, false, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<2, true, false, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<2, true, false, true, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<2, true, true, false, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<2, true, true, false, true>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<2, true, true, true, false>(t,r,mats,profile);
		else if( mode == 2 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<2, true, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<3, false, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<3, false, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<3, false, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<3, false, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<3, false, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<3, false, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<3, false, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<3, false, true, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<3, true, false, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<3, true, false, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<3, true, false, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<3, true, false, true, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<3, true, true, false, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<3, true, true, false, true>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<3, true, true, true, false>(t,r,mats,profile);
		else if( mode == 3 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<3, true, true, true, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<4, false, false, false, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<4, false, false, false, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<4, false, false, true, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<4, false, false, true, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<4, false, true, false, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<4, false, true, false, true>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<4, false, true, true, false>(t,r,mats,profile);
		else if( mode == 4 && !intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<4, false, true, true, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<4, true, false, false, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<4, true, false, false, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<4, true, false, true, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && !intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<4, true, false, true, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && !intv[2] && !privatized)
			mttkrp_combined_lb_5<4, true, true, false, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && !intv[2] && privatized)
			mttkrp_combined_lb_5<4, true, true, false, true>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && intv[2] && !privatized)
			mttkrp_combined_lb_5<4, true, true, true, false>(t,r,mats,profile);
		else if( mode == 4 && intv[0] && intv[1] && intv[2] && privatized)
			mttkrp_combined_lb_5<4, true, true, true, true>(t,r,mats,profile);
	}

	return 0;
}

int b_thread_start(csf* t)
{
	int nmode = t->nmode;
	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	//printf("%d\n",num_th);

	idx_t** bth = new idx_t*[num_th+1];

	idx_t* res = new idx_t[nmode*(num_th+1)];

	for(int i = 0 ; i<num_th + 1; i++)
	{
		bth[i] = res + i*nmode;
	}

	for(int i = 0 ; i < nmode ; i++)
	{	// Each thread is going to go throung nnz/num_th nnz's for even work distribution
		for(int th = 0; th < num_th + 1; th++)
		{
			bth[th][i] = (t->fiber_count[i] * th)/num_th;
		}
	}

	int method = 1;

	if(method == 1)
	{
		#pragma omp parallel for
		for(int th = 1; th < num_th ; th++)
		{	// Search and put into the correct positions
			for(int d = nmode - 1; d > 0 ; d--)
			{	// Do binary search for faster speeds
				while(t->ptr[d-1][bth[th][d-1]] > bth[th][d])
					bth[th][d-1] --;
				
				while(t->ptr[d-1][bth[th][d-1]+1] <= bth[th][d])
					bth[th][d-1] ++;
			}
		}		
	}
	else if (method == 2)
	{
		idx_t** count_tree = new idx_t*[nmode-1];
		for(int i = 0 ; i<nmode-1 ; i++)
		{
			count_tree[i] = new idx_t[ t->fiber_count[i] ];
			//memset(count_tree[i],0,sizeof(int)*t->fiber_count[i]);
		}
		for(int mode = nmode - 2; mode >= 0 ; mode--)
		{
			#ifdef OMP
			#pragma omp parallel for
			#endif
			for(idx_t i = 0 ; i<t->fiber_count[mode] ; i++)
			{
				if(mode == nmode -2)
					count_tree[mode][i] = t->ptr[mode][i+1] - t->ptr[mode][i] + 1; // +1 for self
				else
				{
					count_tree[mode][i] = 1;
					for(int j = t->ptr[mode][i] ; j < t->ptr[mode][i+1] ; j++)
					{
						count_tree[mode][i] += count_tree[mode+1][j];
					}

				}					
			}
		} 
		// Found work per node

		for(int mode = 0; mode < nmode -1 ; mode++)
		{
			// Prefix sum, can be optimized later
			idx_t total_work = 0;
			for(idx_t i = 0 ; i < t->fiber_count[mode] ; i++)
			{
				total_work += count_tree[mode][i];
			}
			
			// Distribute the work between threads
			int th = 1;
			idx_t work = 0;
			for(idx_t i = 0 ; i < t->fiber_count[mode] ; i++)
			{
				work += count_tree[mode][i];
				while (work > (total_work * th)/num_th)
				{
					bth[th][mode] = i;
					th += 1;
				}
			}
		}

		for(int th = 0; th < num_th ; th ++)
		{
			bth[th][nmode-1] = t->ptr[nmode-2][bth[th][nmode-2]];
		}

	}
	
	/*
	for(int i = 0 ; i < nmode ; i++)
	{	// Each thread is going to go throung nnz/num_th nnz's for even work distribution
		for(int th = 0; th < num_th + 1; th++)
		{
			printf("%lld ",bth[th][i]);
		}
		printf("\n");
	}
	*/

	t->b_thread_start = bth;		
	return 0;
}


double atomic_thresh(int r,mutex_array* mutex)
{
	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif
	double thresh = 1;
	if (num_th == 1)
		return thresh;
	
	idx_t small_array_size = 10000;
	idx_t max_search_limit = 128;
	idx_t min_coef = 1;
	idx_t max_coef = max_search_limit;
	idx_t coef = min_coef;
	idx_t big_array_size = coef*small_array_size;
	bool atomic_better = true;
	TYPE* res_array = new TYPE[small_array_size*r];
	TYPE** priv_array = new TYPE*[num_th];
	for(int i=0; i<num_th ; i++)
	{
		priv_array[i] = new TYPE[small_array_size*r];
	}
	
	idx_t* rand_array  = new idx_t[small_array_size*max_search_limit];
	srand(0);
	for(idx_t i=0; i<small_array_size*max_search_limit ; i++)
	{
		rand_array[i] = rand()%small_array_size;
	}
	while(atomic_better &&  max_coef > min_coef + 1)
	{
		double atomic_time, privatized_time;
		{
			auto start = std::chrono::high_resolution_clock::now();
			// Do atomic computation
			#pragma omp parallel for
			for(idx_t i=0 ; i< big_array_size ; i++)
			{
				idx_t row_id = rand_array[i];
				#ifdef OMP
				mutex_set_lock(mutex,row_id);
				#endif
				#pragma omp simd
				for(int y = 0 ; y<r ; y++)
				{	
					res_array[row_id*r + y] += i*y;
				}
				#ifdef OMP
				mutex_unset_lock(mutex,row_id);
				#endif
			}
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			atomic_time = diff.count();
		}	
		{
			auto start = std::chrono::high_resolution_clock::now();
			// Do privatized computation
			#pragma omp parallel
			{
				#ifdef OMP
				int th = omp_get_thread_num();
				#else
				int th = 0;
				#endif
				#pragma omp for
				for(idx_t i=0 ; i< big_array_size ; i++)
				{
					idx_t row_id = rand_array[i];
					#pragma omp simd
					for(int y = 0 ; y<r ; y++)
					{	
						priv_array[th][row_id*r + y] += i*y;
					}
				}
			}
			// reduction step
			#pragma omp parallel for
			for(idx_t i = 0 ; i<small_array_size ; i++)
			{
				for(int th = 0; th<num_th ; th++)
				{
					for(int y=0; y<r ;y++)
						res_array[i*r + y] += priv_array[th][i*r + y];
				}
			}
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			privatized_time = diff.count();
		}

		if(privatized_time < atomic_time)
		{
			max_coef = coef;
		}
		else
		{
			min_coef = coef;
		}
		
		if(max_coef > min_coef + 1)
		{
			coef = (max_coef + min_coef ) / 2;
			big_array_size = coef*small_array_size;
		}		
	}

	thresh = 1/((double) coef);
	

	return thresh;
}

int reduce_mode_0(csf* t, matrix* mat)
{
	#pragma omp barrier

	idx_t** thread_start = t->b_thread_start;
	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif
	// Handle boundaries first
	for(int th=num_th -1 ; th>0 ; th--)
	{
		idx_t loc = t->ind[0][ thread_start[th][0] ];
		TYPE* target = mat->val + ( loc + th - 1)*(mat->dim2);
		TYPE* source = mat->val + (loc + th)* (mat->dim2);
		for(int y = 0 ; y<mat->dim2 ; y++)
		{
			target[y] += source[y];
			source[y] = 0;
		}
	}
	
	for(int th = 1; th < num_th ; th++)
	{
		idx_t a = thread_start[th+1][0]+1;
		idx_t b = t->fiber_count[0];
		idx_t end = MIN(a,b);
		idx_t start = thread_start[th][0]+1;
		for(idx_t i=start; i<end; i ++)
		{
			TYPE* target = mat->val + (t->ind[0][i])*mat->dim2;
			TYPE* source = mat->val + (t->ind[0][i] + th)*mat->dim2;
			for(int y = 0 ; y<mat->dim2 ; y++)
			{
				target[y] = source[y];
				source[y] = 0;
			}
		}
	}
	return 0;
}


int reduce_socket(csf* t, int r, matrix* mat)
{
	#pragma omp parallel for	
	for(int i = 0 ; i < mat->dim1 ; i ++)
	{
		TYPE* out = mat->val + (mat->dim2)*i;
		for(int socket = 1 ; socket < t->num_sockets ; socket++)
		{
			TYPE* in = t->private_mats[socket*t->cps]->val + (mat->dim2)*i;
			for(int j = 0 ; j < r ; j++)
			{
				out[j] += in[j];
			}
		}
	}

	return 0;
}

#endif