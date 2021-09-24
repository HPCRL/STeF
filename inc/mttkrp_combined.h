#ifndef MTTKRP_COMBINED
#define MTTKRP_COMBINED

#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"
#include "../inc/mutex.h"
#include "../inc/mttkrp_hardwired.h"


template <int mode, bool intv1, bool privatized>
int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile )
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
	
	if (mode > 0 && privatized)
	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	
	mutex_array* mutex = t->mutex;

	//set_matrix(*mats[mode],0);
	memset(mats[mode]->val,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
	
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
						else if (privatized)
						{                            
							TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								matval[y] += tval * pr[y];	// saxpy step			
							}
						}
						else
						{
							const int row_id = t->ind[2][i2];
							#ifdef OMP
							mutex_set_lock(mutex,row_id);
							#endif
							TYPE* matval = mats[2]->val + ((mats[2]) -> dim2) * row_id;
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								matval[y] += tval * pr[y];	// saxpy step			
							}
							#ifdef OMP
							mutex_unset_lock(mutex,row_id);
							#endif
						}
					}	
				}			
				
				if(mode == 0)
				{
					//TYPE* matval = mv1;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						mv1[y] += pr[y] * mv2[y]; // dot TTV
					}
				}
				else if (mode == 1)
				{
					if(privatized)
					{
						TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr[y] * mv1[y]; // saxpy
						}
					}
					else
					{
						const int row_id = t->ind[1][i1];
						#ifdef OMP
						mutex_set_lock(mutex,row_id);
						#endif

						TYPE* matval = mats[1]->val + ((mats[1]) -> dim2) * row_id;
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr[y] * mv1[y]; // saxpy
						}

						#ifdef OMP
						mutex_unset_lock(mutex,row_id);
						#endif
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
	
	
	if (mode > 0 && privatized)
	{
		reduce(t,r,mats[mode]);
	}
	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}

int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile , int mode, bool intv1 );

template <int mode, bool intv1, bool intv2>
int mttkrp_combined_4(csf* t, int r, matrix** mats, int profile );

template <int mode, bool intv1, bool intv2, bool intv3>
int mttkrp_combined_5(csf* t, int r, matrix** mats, int profile );

int mttkrp_combined(csf* t, int r, matrix** mats, int profile , int mode, bool privatized, bool* intv );

#endif