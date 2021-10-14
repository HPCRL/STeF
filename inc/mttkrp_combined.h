#ifndef MTTKRP_COMBINED
#define MTTKRP_COMBINED

#include "../inc/matrix.h"
#include <stdlib.h> 
#include "../inc/tensor.h"
#include "../inc/util.h"
#include "../inc/mutex.h"
#include "../inc/mttkrp_hardwired.h"

int mttkrp_combined(csf* t, int r, matrix** mats, int profile , int mode, bool privatized, bool* intv );
int mttkrp_combined_lb(csf* t, int r, matrix** mats, int profile , int mode, bool privatized, bool* intv );
int b_thread_start(csf* t);
double atomic_thresh(int r,mutex_array* mutex);
int reduce_mode_0(csf* t, matrix* mat);

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
		TYPE* pr = partial_results;
		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			TYPE* mv1 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				
				TYPE* mv2 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if (intv1)
				{
					pr = t->intval[1] + i1*r ;
					//memset(pr,0,sizeof(TYPE)*r);
				}

				if(mode == 0 || !(intv1 && mode<2))
					memset(pr,0,sizeof(TYPE)*r);

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
	if(profile == mode)
	{
		LIKWID_MARKER_CLOSE;
	}
	rem(partial_results_all);	
	return 0;
	
}

int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile , int mode, bool intv1 );

template <int mode, bool intv1, bool intv2, bool privatized>
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
//	TYPE* dot_partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	{
		partial_results_all[i] = 0;
	}
		
	
	if (mode > 0 && privatized)
	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	
	mutex_array* mutex = t->mutex;

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
		// TYPE* dot_partial_results = dot_partial_results_all + th*partial_results_size;

		TYPE* pr0 = partial_results;
		TYPE* pr1 = partial_results + r;


		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			TYPE* mv0 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* mv1 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if(intv1 && mode < 2)
					pr0 = t->intval[1] + i1*r ; // Set the location for intermediate value for T(i_1,i_2)
				if(mode == 0 || !(intv1 && mode<2))
					memset(pr0,0,sizeof(TYPE)*r); // Reset the previous values otherwise

				if(mode >= 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr0[y] += mv0[y] * mv1[y];	// KrP step			
					}
				}

				if( mode == 0 || mode > 1 || !intv1 )
				{					
					for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
					{	
						TYPE* mv2 = mats[2]->val + ((mats[2]) -> dim2) * t->ind[2][i2];						
						if(intv2 && mode < 3)
							pr1 = t->intval[2] + i2*r; // Set the location for intermediate value for T(i_1,i_2,i_3)
						if(mode == 0 || !(intv2 && mode<3))
							memset(pr1,0,sizeof(TYPE)*r); // Reset the previous values otherwise									

						if(mode > 2)
						{
							for(int y=0 ; y<r ; y++)
							{
								pr1[y] += pr0[y] * mv2[y] ;	// KrP step			
							}
						}

						if( mode == 0 || mode > 2 || !intv2 )
						{					
							for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
							{	
								TYPE* mv3 = mats[3]->val + ((mats[3] -> dim2) * t->ind[3][i3]);
								TYPE tval = t->val[i3];													
							
								if(mode < 3)
								{				
									// Do dot product TTM for mode 3					
									TYPE* matval = mv3;
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr1[y] += tval * matval[y];	// dot TTM step			
									}
								}
								else if (privatized)
								{                            
									TYPE* matval = t->private_mats[th]->val + ((mats[3]) -> dim2) * t->ind[3][i3];
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										matval[y] += tval * pr1[y];	// saxpy step			
									}
								}
								else
								{
									TYPE* matval = mv3;
									idx_t row_id = t->ind[3][i3] ; 
									#ifdef OMP
									mutex_set_lock(mutex,row_id);
									#endif									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										matval[y] += tval * pr1[y];	// saxpy step			
									}
									#ifdef OMP
									mutex_unset_lock(mutex,row_id);
									#endif
								}
							}
						} // End of loop for traversing mode 3


						if(mode < 2)
						{	// Do dot product TTV for mode 2							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr0[y] += pr1[y] * mv2[y];	// dot TTV step			
							}
						}
						else if(mode == 2)
						{	// MTTKRP for mode 2					
							if (privatized)
							{
								TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
							}
							else
							{
								idx_t row_id = t->ind[2][i2] ; 
								TYPE* matval = mv2;
								#ifdef OMP
								mutex_set_lock(mutex,row_id);
								#endif
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
								#ifdef OMP
								mutex_unset_lock(mutex,row_id);
								#endif
							}								
						}

					}	
				} // End of the loop for traversing mode 2	
				
				if(mode < 1)
				{
					//TYPE* matval = mv1;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						mv0[y] += pr0[y] * mv1[y]; // dot TTV to compute mode 0
					}
				}
				else if (mode == 1)
				{	// MTTKRP for mode 1
					if(privatized)
					{
						TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
					}
					else
					{
						TYPE* matval = mv1;
						idx_t row_id = t->ind[1][i1] ; 
						#ifdef OMP
						mutex_set_lock(mutex,row_id);
						#endif
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
						#ifdef OMP
						mutex_unset_lock(mutex,row_id);
						#endif
					}
					
				}
			} // End of loop for traversing mode 1
		} // End of loop for traversing mode 0
		


		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	} // End of parallel region
	
	
	if (mode > 0 && privatized)
	{
		reduce(t,r,mats[mode]);
	}

	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}

template <int mode, bool intv1, bool intv2, bool intv3, bool privatized>
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
	{
		partial_results_all[i] = 0;
	}
		
	if (mode > 0 && privatized)
	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	
	mutex_array* mutex = t->mutex;

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
		TYPE* pr0 = partial_results;
		TYPE* pr1 = partial_results + r;
		TYPE* pr2 = partial_results + 2*r;


		for(idx_t i0 = thread_start[th] ; i0 < thread_start[th+1] ; i0++)
		{
			TYPE* mv0 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			for(idx_t i1 = t->ptr[0][i0] ; i1< t->ptr[0][i0+1]; i1++)
			{
				TYPE* mv1 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if(intv1 && mode < 2)
					pr0 = t->intval[1] + i1*r ; // Set the location for intermediate value for T(i_1,i_2)
				if(mode == 0 || !(intv1 && mode<2))
					memset(pr0,0,sizeof(TYPE)*r); // Reset the previous values otherwise


				if(mode >= 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr0[y] += mv0[y] * mv1[y];	// KrP step			
					}
				}

				if( mode == 0 || mode > 1 || !intv1 )
				{					
					for(idx_t i2 = t->ptr[1][i1] ; i2< t->ptr[1][i1+1]; i2++)
					{	
						TYPE* mv2 = mats[2]->val + ((mats[2]) -> dim2) * t->ind[2][i2];						
						if(intv2 && mode < 3)
							pr1 = t->intval[2] + i2*r; // Set the location for intermediate value for T(i_1,i_2,i_3)
						if(mode == 0 || !(intv2 && mode<3))
							memset(pr1,0,sizeof(TYPE)*r); // Reset the previous values otherwise						

						if(mode > 2)
						{
							for(int y=0 ; y<r ; y++)
							{
								pr1[y] += pr0[y] * mv2[y] ;	// KrP step			
							}
						}

						if( mode == 0 || mode > 2 || !intv2 )
						{					
							for(idx_t i3 = t->ptr[2][i2] ; i3< t->ptr[2][i2+1]; i3++)
							{	
								TYPE* mv3 = mats[3]->val + ((mats[3] -> dim2) * t->ind[3][i3]);
								if(intv3 && mode < 4)
									pr2 = t->intval[3] + i3*r; // Set the location for intermediate value for T(i_1,i_2,i_3,i_4)
								if(mode == 0 || !(intv3 && mode<4))
									memset(pr2,0,sizeof(TYPE)*r); // Reset the previous values otherwise	
								//TYPE tval = t->val[i3];													
							
								if(mode > 3)
								{
									for(int y=0 ; y<r ; y++)
									{
										pr2[y] += pr1[y] * mv3[y] ;	// KrP step			
									}
								}

								if( mode == 0 || mode > 3 || !intv3 )
								{					
									for(idx_t i4 = t->ptr[3][i3] ; i4< t->ptr[3][i3+1]; i4++)
									{	
										TYPE* mv4 = mats[4]->val + ((mats[4] -> dim2) * t->ind[4][i4]);
										TYPE tval = t->val[i4];	

										if(mode < 4)
										{				
											// Do dot product TTM for mode 3					
											TYPE* matval = mv4;
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												pr2[y] += tval * matval[y];	// dot TTM step			
											}
										}
										else if (mode == 4)
										{								
											if(privatized)
											{                            
												TYPE* matval = t->private_mats[th]->val + ((mats[4]) -> dim2) * t->ind[4][i4];
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													matval[y] += tval * pr2[y];	// saxpy step			
												}
											}
											else
											{
												TYPE* matval = mv4;
												idx_t row_id = t->ind[4][i4] ; 
												#ifdef OMP
												mutex_set_lock(mutex,row_id);
												#endif									
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													matval[y] += tval * pr2[y];	// saxpy step			
												}
												#ifdef OMP
												mutex_unset_lock(mutex,row_id);
												#endif
											}
										}

									}
								}

								if(mode < 3)
								{				
									// Do dot product TTM for mode 3					
									TYPE* matval = mv3;
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr1[y] += pr2[y] * matval[y];	// dot TTV step			
									}
								}
								else if (mode == 3)
								{								
									if(privatized)
									{                            
										TYPE* matval = t->private_mats[th]->val + ((mats[3]) -> dim2) * t->ind[3][i3];
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											matval[y] += pr1[y] * pr2[y];	// saxpy step			
										}
									}
									else
									{
										TYPE* matval = mv3;
										idx_t row_id = t->ind[3][i3] ; 
										#ifdef OMP
										mutex_set_lock(mutex,row_id);
										#endif									
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											matval[y] += pr1[y] * pr2[y];	// saxpy step			
										}
										#ifdef OMP
										mutex_unset_lock(mutex,row_id);
										#endif
									}
								}
							}
						} // End of loop for traversing mode 3


						if(mode < 2)
						{	// Do dot product TTV for mode 2							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr0[y] += pr1[y] * mv2[y];	// dot TTV step			
							}
						}
						else if(mode == 2)
						{	// MTTKRP for mode 2					
							if (privatized)
							{
								TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
							}
							else
							{
								idx_t row_id = t->ind[2][i2] ; 
								TYPE* matval = mv2;
								#ifdef OMP
								mutex_set_lock(mutex,row_id);
								#endif
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
								#ifdef OMP
								mutex_unset_lock(mutex,row_id);
								#endif
							}								
						}

					}	
				} // End of the loop for traversing mode 2	
				
				if(mode < 1)
				{
					//TYPE* matval = mv1;
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						mv0[y] += pr0[y] * mv1[y]; // dot TTV to compute mode 0
					}
				}
				else if (mode == 1)
				{	// MTTKRP for mode 1
					if(privatized)
					{
						TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
					}
					else
					{
						TYPE* matval = mv1;
						idx_t row_id = t->ind[1][i1] ; 
						#ifdef OMP
						mutex_set_lock(mutex,row_id);
						#endif
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
						#ifdef OMP
						mutex_unset_lock(mutex,row_id);
						#endif
					}
					
				}
			} // End of loop for traversing mode 1
		} // End of loop for traversing mode 0
		


		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	} // End of parallel region
	
	
	if (mode > 0 && privatized)
	{
		reduce(t,r,mats[mode]);
	}

	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	}

template <int mode, bool intv1, bool privatized>
int mttkrp_combined_lb_3(csf* t, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; 
	idx_t** thread_start = t->b_thread_start;
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

		TYPE* pr = partial_results;

		idx_t i0_end = MIN(thread_start[th+1][0] + 1 , t->fiber_count[0]);
		for(idx_t i0 = thread_start[th][0] ; i0 < i0_end ; i0++)
		{
			TYPE* mv1 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];
			idx_t i1_start = MAX(t->ptr[0][i0],thread_start[th][1]);
			idx_t i1_end = MIN(t->ptr[0][i0+1],thread_start[th+1][1]+1);
			for(idx_t i1 = i1_start ; i1< i1_end; i1++)
			{
				TYPE* mv2 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];				
				
				
				if(intv1 && mode < 2)
				{
					pr = t->intval[1] + (i1+th)*r;
				}
				if(mode == 0 || !(intv1 && mode<2))
					memset(pr,0,sizeof(TYPE)*r);
					

				if(mode == 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr[y] += mv1[y] * mv2[y];	// KrP step			
					}
				}

				if( mode != 1 || !intv1 )
				{
					idx_t i2_start = MAX(t->ptr[1][i1],thread_start[th][2]);
					idx_t i2_end = MIN(t->ptr[1][i1+1],thread_start[th+1][2]);
					for(idx_t i2 = i2_start ; i2< i2_end; i2++)
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
					TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * (t->ind[0][i0]+th);
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						matval[y] += pr[y] * mv2[y]; // dot TTV
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
	if(mode == 0)
	{
		reduce_mode_0(t,mats[0]);
	}

	if (mode > 0 && privatized)
	{
		reduce(t,r,mats[mode]);
	}
	if(profile == mode)
	{
		LIKWID_MARKER_CLOSE;
	}
	rem(partial_results_all);	
	return 0;
	
}

template <int mode, bool intv1, bool intv2, bool privatized>
int mttkrp_combined_lb_4(csf* t, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; 
	idx_t** thread_start = t->b_thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
//	TYPE* dot_partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	{
		partial_results_all[i] = 0;
	}
		
	
	if (mode > 0 && privatized)	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	
	mutex_array* mutex = t->mutex;

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
		// TYPE* dot_partial_results = dot_partial_results_all + th*partial_results_size;

		TYPE* pr0 = partial_results;
		TYPE* pr1 = partial_results + r;

		idx_t i0_end = MIN(thread_start[th+1][0]+1,t->fiber_count[0]);
		for(idx_t i0 = thread_start[th][0] ; i0 < i0_end ; i0++)
		{
			TYPE* mv0 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];

			idx_t i1_start = MAX(t->ptr[0][i0],thread_start[th][1]);
			idx_t i1_end = MIN(t->ptr[0][i0+1],thread_start[th+1][1]+1);
			for(idx_t i1 = i1_start ; i1< i1_end; i1++)
			{
				TYPE* mv1 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if(intv1 && mode < 2)
					pr0 = t->intval[1] + (th + i1)*r ; // Set the location for intermediate value for T(i_1,i_2)
				if(mode == 0 || !(intv1 && mode<2))
					memset(pr0,0,sizeof(TYPE)*r); // Reset the previous values otherwise


				if(mode >= 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr0[y] += mv0[y] * mv1[y];	// KrP step			
					}
				}

				if( mode == 0 || mode > 1 || !intv1 )
				{		
					idx_t i2_start = MAX(t->ptr[1][i1],thread_start[th][2]);
					idx_t i2_end = MIN(t->ptr[1][i1+1],thread_start[th+1][2]+1);			
					for(idx_t i2 = i2_start ; i2< i2_end; i2++)
					{	
						TYPE* mv2 = mats[2]->val + ((mats[2]) -> dim2) * t->ind[2][i2];						
						if(intv2 && mode < 3)
							pr1 = t->intval[2] + (th + i2)*r; // Set the location for intermediate value for T(i_1,i_2,i_3)
						if(mode == 0 || !(intv2 && mode<3))
							memset(pr1,0,sizeof(TYPE)*r); // Reset the previous values otherwise							

						if(mode > 2)
						{
							for(int y=0 ; y<r ; y++)
							{
								pr1[y] += pr0[y] * mv2[y] ;	// KrP step			
							}
						}

						if( mode == 0 || mode > 2 || !intv2 )
						{					
							idx_t i3_start = MAX(t->ptr[2][i2],thread_start[th][3]);
							idx_t i3_end = MIN(t->ptr[2][i2+1],thread_start[th+1][3]);			
							for(idx_t i3 = i3_start ; i3< i3_end; i3++)
							{	
								TYPE* mv3 = mats[3]->val + ((mats[3] -> dim2) * t->ind[3][i3]);
								TYPE tval = t->val[i3];													
							
								if(mode < 3)
								{				
									// Do dot product TTM for mode 3					
									TYPE* matval = mv3;
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr1[y] += tval * matval[y];	// dot TTM step			
									}
								}
								else if (privatized)
								{                            
									TYPE* matval = t->private_mats[th]->val + ((mats[3]) -> dim2) * t->ind[3][i3];
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										matval[y] += tval * pr1[y];	// saxpy step			
									}
								}
								else
								{
									TYPE* matval = mv3;
									idx_t row_id = t->ind[3][i3] ; 
									#ifdef OMP
									mutex_set_lock(mutex,row_id);
									#endif									
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										matval[y] += tval * pr1[y];	// saxpy step			
									}
									#ifdef OMP
									mutex_unset_lock(mutex,row_id);
									#endif
								}
							}
						} // End of loop for traversing mode 3


						if(mode < 2)
						{	// Do dot product TTV for mode 2							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr0[y] += pr1[y] * mv2[y];	// dot TTV step			
							}
						}
						else if(mode == 2)
						{	// MTTKRP for mode 2					
							if (privatized)
							{
								TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
							}
							else
							{
								idx_t row_id = t->ind[2][i2] ; 
								TYPE* matval = mv2;
								#ifdef OMP
								mutex_set_lock(mutex,row_id);
								#endif
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
								#ifdef OMP
								mutex_unset_lock(mutex,row_id);
								#endif
							}								
						}

					}	
				} // End of the loop for traversing mode 2	
				
				if(mode < 1)
				{
					TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * (t->ind[0][i0]+th);
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						matval[y] += pr0[y] * mv1[y]; // dot TTV
					}	
				}
				else if (mode == 1)
				{	// MTTKRP for mode 1
					if(privatized)
					{
						TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // dot TTV
						}
					}
					else
					{
						TYPE* matval = mv1;
						idx_t row_id = t->ind[1][i1] ; 
						#ifdef OMP
						mutex_set_lock(mutex,row_id);
						#endif
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
						#ifdef OMP
						mutex_unset_lock(mutex,row_id);
						#endif
					}
					
				}
			} // End of loop for traversing mode 1
		} // End of loop for traversing mode 0
		


		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	} // End of parallel region

	if(mode == 0)
	{
		reduce_mode_0(t,mats[0]);
	}

	if (mode > 0 && privatized)	
	{
		reduce(t,r,mats[mode]);
	}

	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
	
}


template <int mode, bool intv1, bool intv2, bool intv3, bool privatized>
int mttkrp_combined_lb_5(csf* t, int r, matrix** mats, int profile )
{
	int nmode = t->nmode;
	int num_th = 1;
	int partial_results_size = nmode*r+PAD; 
	idx_t** thread_start = t->b_thread_start;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif
	
	if(VERBOSE >= VERBOSE_DEBUG) printf("num ths %d\n", num_th);
	
	TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));
	
	for(int i = 0 ; i < num_th*partial_results_size ; i++)
	{
		partial_results_all[i] = 0;
	}
		
	if (mode > 0 && privatized)	{
		for ( int i=0; i<num_th; i++)
		{
			memset(t->private_mats[i]->val ,0,sizeof(TYPE)*mats[mode]->dim1 * mats[mode]->dim2);
		}
	}
	
	mutex_array* mutex = t->mutex;

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
		TYPE* pr0 = partial_results;
		TYPE* pr1 = partial_results + r;
		TYPE* pr2 = partial_results + 2*r;


		idx_t i0_end = MIN(thread_start[th+1][0]+1,t->fiber_count[0]);
		for(idx_t i0 = thread_start[th][0] ; i0 < i0_end ; i0++)
		{
			TYPE* mv0 = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0];

			idx_t i1_start = MAX(t->ptr[0][i0],thread_start[th][1]);
			idx_t i1_end = MIN(t->ptr[0][i0+1],thread_start[th+1][1]+1);
			for(idx_t i1 = i1_start ; i1< i1_end; i1++)
			{
				TYPE* mv1 = mats[1]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
				if(intv1 && mode < 2)
					pr0 = t->intval[1] + (th + i1)*r ; // Set the location for intermediate value for T(i_1,i_2)
				if(mode == 0 || !(intv1 && mode<2))
					memset(pr0,0,sizeof(TYPE)*r); // Reset the previous values otherwise


				if(mode >= 2)
				{
					for(int y=0 ; y<r ; y++)
					{
						pr0[y] += mv0[y] * mv1[y];	// KrP step			
					}
				}

				if( mode == 0 || mode > 1 || !intv1 )
				{					
					idx_t i2_start = MAX(t->ptr[1][i1],thread_start[th][2]);
					idx_t i2_end = MIN(t->ptr[1][i1+1],thread_start[th+1][2]+1);
					for(idx_t i2 = i2_start ; i2< i2_end; i2++)
					{	
						TYPE* mv2 = mats[2]->val + ((mats[2]) -> dim2) * t->ind[2][i2];						
						if(intv2 && mode < 3)
							pr1 = t->intval[2] + (th + i2)*r; // Set the location for intermediate value for T(i_1,i_2,i_3)
						if(mode == 0 || !(intv2 && mode<3))
							memset(pr1,0,sizeof(TYPE)*r); // Reset the previous values otherwise							

						if(mode > 2)
						{
							for(int y=0 ; y<r ; y++)
							{
								pr1[y] += pr0[y] * mv2[y] ;	// KrP step			
							}
						}

						if( mode == 0 || mode > 2 || !intv2 )
						{					
							idx_t i3_start = MAX(t->ptr[2][i2],thread_start[th][3]);
							idx_t i3_end = MIN(t->ptr[2][i2+1],thread_start[th+1][3]+1);
							for(idx_t i3 = i3_start ; i3< i3_end; i3++)
							{	
								TYPE* mv3 = mats[3]->val + ((mats[3] -> dim2) * t->ind[3][i3]);
								if(intv3 && mode < 4)
									pr2 = t->intval[3] + (th + i3)*r; // Set the location for intermediate value for T(i_1,i_2,i_3,i_4)
								if(mode == 0 || !(intv3 && mode<4))
									memset(pr2,0,sizeof(TYPE)*r); // Reset the previous values otherwise	
								//TYPE tval = t->val[i3];													
							
								if(mode > 3)
								{
									for(int y=0 ; y<r ; y++)
									{
										pr2[y] += pr1[y] * mv3[y] ;	// KrP step			
									}
								}

								if( mode == 0 || mode > 3 || !intv3 )
								{					
									idx_t i4_start = MAX(t->ptr[3][i3],thread_start[th][4]);
									idx_t i4_end = MIN(t->ptr[3][i3+1],thread_start[th+1][4]);
									for(idx_t i4 = i4_start ; i4< i4_end; i4++)
									{	
										TYPE* mv4 = mats[4]->val + ((mats[4] -> dim2) * t->ind[4][i4]);
										TYPE tval = t->val[i4];	

										if(mode < 4)
										{				
											// Do dot product TTM for mode 3					
											TYPE* matval = mv4;
											#pragma omp simd
											for(int y=0 ; y<r ; y++)
											{
												pr2[y] += tval * matval[y];	// dot TTM step			
											}
										}
										else if (mode == 4)
										{								
											if(privatized)
											{                            
												TYPE* matval = t->private_mats[th]->val + ((mats[4]) -> dim2) * t->ind[4][i4];
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													matval[y] += tval * pr2[y];	// saxpy step			
												}
											}
											else
											{
												TYPE* matval = mv4;
												idx_t row_id = t->ind[4][i4] ; 
												#ifdef OMP
												mutex_set_lock(mutex,row_id);
												#endif									
												#pragma omp simd
												for(int y=0 ; y<r ; y++)
												{
													matval[y] += tval * pr2[y];	// saxpy step			
												}
												#ifdef OMP
												mutex_unset_lock(mutex,row_id);
												#endif
											}
										}

									}
								}

								if(mode < 3)
								{				
									// Do dot product TTM for mode 3					
									TYPE* matval = mv3;
									#pragma omp simd
									for(int y=0 ; y<r ; y++)
									{
										pr1[y] += pr2[y] * matval[y];	// dot TTV step			
									}
								}
								else if (mode == 3)
								{								
									if(privatized)
									{                            
										TYPE* matval = t->private_mats[th]->val + ((mats[3]) -> dim2) * t->ind[3][i3];
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											matval[y] += pr1[y] * pr2[y];	// saxpy step			
										}
									}
									else
									{
										TYPE* matval = mv3;
										idx_t row_id = t->ind[3][i3] ; 
										#ifdef OMP
										mutex_set_lock(mutex,row_id);
										#endif									
										#pragma omp simd
										for(int y=0 ; y<r ; y++)
										{
											matval[y] += pr1[y] * pr2[y];	// saxpy step			
										}
										#ifdef OMP
										mutex_unset_lock(mutex,row_id);
										#endif
									}
								}
							}
						} // End of loop for traversing mode 3


						if(mode < 2)
						{	// Do dot product TTV for mode 2							
							#pragma omp simd
							for(int y=0 ; y<r ; y++)
							{
								pr0[y] += pr1[y] * mv2[y];	// dot TTV step			
							}
						}
						else if(mode == 2)
						{	// MTTKRP for mode 2					
							if (privatized)
							{
								TYPE* matval = t->private_mats[th]->val + ((mats[2]) -> dim2) * t->ind[2][i2];
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
							}
							else
							{
								idx_t row_id = t->ind[2][i2] ; 
								TYPE* matval = mv2;
								#ifdef OMP
								mutex_set_lock(mutex,row_id);
								#endif
								#pragma omp simd
								for(int y=0 ; y<r ; y++)
								{
									matval[y] += pr0[y] * pr1[y];	// saxpy step			
								}
								#ifdef OMP
								mutex_unset_lock(mutex,row_id);
								#endif
							}								
						}

					}	
				} // End of the loop for traversing mode 2	
				
				if(mode < 1)
				{
					TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * (t->ind[0][i0]+th);
					#pragma omp simd
					for(int y=0 ; y<r ; y++)
					{
						matval[y] += pr0[y] * mv1[y]; // dot TTV
					}	
				}
				else if (mode == 1)
				{	// MTTKRP for mode 1
					if(privatized)
					{
						TYPE* matval = t->private_mats[th]->val + ((mats[1]) -> dim2) * t->ind[1][i1];
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
					}
					else
					{
						TYPE* matval = mv1;
						idx_t row_id = t->ind[1][i1] ; 
						#ifdef OMP
						mutex_set_lock(mutex,row_id);
						#endif
						#pragma omp simd
						for(int y=0 ; y<r ; y++)
						{
							matval[y] += pr0[y] * mv0[y]; // saxpy
						}
						#ifdef OMP
						mutex_unset_lock(mutex,row_id);
						#endif
					}
					
				}
			} // End of loop for traversing mode 1
		} // End of loop for traversing mode 0
		


		auto time_end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_diff = time_end-time_start;
		
		
		if(VERBOSE >= VERBOSE_DEBUG) printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		
		if(profile == mode)
		{
			LIKWID_MARKER_STOP("Compute");
		}
	} // End of parallel region
	
	if(mode == 0)
	{
		reduce_mode_0(t,mats[0]);
	}
	
	if (mode > 0 && privatized)	
	{
		reduce(t,r,mats[mode]);
	}

	LIKWID_MARKER_CLOSE;
	rem(partial_results_all);	
	return 0;
}


#endif