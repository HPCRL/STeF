#ifndef MTTKRP_HARDWIRED_CPP
#define MTTKRP_HARDWIRED_CPP
#include "../inc/mttkrp_hardwired.h"

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
	}
	rem(partial_results_all);	
	return 0;

}

#endif