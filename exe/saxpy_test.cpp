#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/util.h"
#include <iostream>

using namespace std;

struct Reverse_Indices
{
	idx_t* krp;
	idx_t* inds;
	TYPE** vals;
	idx_t* ptr;
	idx_t tile_size
};


Reverse_Indices* find_reverse_indices(csf* t, idx_t tile )
{
	Reverse_Indices* res = new Reverse_Indices[1];
        Reverse_Indices ri = *res;
        int mode = t->nmode-1;
        int mlen = t->mlen[mode];
        ri.ptr = new idx_t[mlen+1];
        idx_t fibcnt = t->fiber_count[mode];
        idx_t* cnt = new idx_t[mlen];
        ri.vals = new TYPE*[fibcnt];
        ri.krp = new idx_t[fibcnt];
	ri.inds = new idx_t[fibcnt];
	ri.tile_size = tile;

        memset(cnt,0,sizeof(idx_t)*(mlen));
        memset(ri.ptr,0,sizeof(idx_t)*(mlen+1));

	idx_t num_tiles = (t->fiber_count[mode-1] - 1) / ri.tile_size + 1;

	
        for(idx_t i=0 ; i< fibcnt; i++)
        {
                ri.ptr[t->ind[mode][i]+1] ++;
	}

        for(idx_t i=0; i<mlen; i++)
        {
                ri.ptr[i+1] += ri.ptr[i];
                cnt[i] = ri.ptr[i];
        }

        cout<<ri.ptr[mlen]<<endl;


        for(idx_t i = 0; i< t->fiber_count[mode-1]; i++)
        {
                for(idx_t j=t->ptr[mode-1][i] ; j< t->ptr[mode-1][i+1]; j++)
                {
                        idx_t row = t->ind[mode][j];
                        idx_t loc = cnt[row];
                        cnt[row]++;
                        ri.krp[loc] = i;
                        ri.vals[loc] = t->val + j;

                }
        }

        *res = ri;
        return res;
}
Reverse_Indices* find_reverse_indices(csf* t)
{
	Reverse_Indices* res = new Reverse_Indices[1];
	Reverse_Indices ri = *res;
	int mode = t->nmode-1;
	int mlen = t->mlen[mode];
	ri.ptr = new idx_t[mlen+1];
	idx_t fibcnt = t->fiber_count[mode];
	idx_t* cnt = new idx_t[mlen];
	ri.vals = new TYPE*[fibcnt];
	ri.krp = new idx_t[fibcnt];
	memset(cnt,0,sizeof(idx_t)*(mlen));
	memset(ri.ptr,0,sizeof(idx_t)*(mlen+1));
	
	for(idx_t i=0 ; i< t-> fiber_count[mode]; i++)
	{
		ri.ptr[t->ind[mode][i]+1] ++;
	}
	
	for(idx_t i=0; i<mlen; i++)
	{
		ri.ptr[i+1] += ri.ptr[i];
		cnt[i] = ri.ptr[i];
	}
	
	cout<<ri.ptr[mlen]<<endl;

	
	for(idx_t i = 0; i< t->fiber_count[mode-1]; i++)
	{
		for(idx_t j=t->ptr[mode-1][i] ; j< t->ptr[mode-1][i+1]; j++)
		{
			idx_t row = t->ind[mode][j];
			idx_t loc = cnt[row];
			cnt[row]++;
			ri.krp[loc] = i;
			ri.vals[loc] = t->val + j;
			
		}
	}
	
	*res = ri;
	return res;
}

int saxpy(csf * t, matrix** mats,int r)
{
	const idx_t nmodes = t->nmode;
	idx_t mode = nmodes-1;
	
	memset(mats[mode]->val, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		#ifdef OMP
		int th = omp_get_thread_num();
		#else
		int th = 0;
		#endif
		idx_t start  = t->thread_start[th];
		idx_t stop = t->thread_start[th+1];
		TYPE* const part_val = new TYPE[nmodes * r];
		idx_t* stack = new idx_t[nmodes];
		idx_t depth = 0;
		for(idx_t s=start ; s<stop ; s++)
		{
			stack[0] = s;
			idx_t sid = t->ind[0][s];
			const TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * sid; 
			
			for(int i=0 ; i< r; i++)
			{
				part_val[i] = matval[i];
			}
			depth ++;
			idx_t pos=s;
			for(int i=0; i<nmodes-2 ; i++)
			{
				pos = t->ptr[i][pos];
			}
			for(int i=1; i<nmodes ; i++)
			{
				stack[i] = t->ptr[i-1][stack[i-1]];
				
			}
			
			
			while( depth > 0 && stack[0] == s)
			{
				for(;depth < mode; depth++)
				{
					idx_t sfid = t->ind[depth][stack[depth]];
					const TYPE* matvali = mats[depth]->val + ((mats[depth]) -> dim2) * sfid; 
					TYPE * const target = part_val + depth*r;
					const TYPE* source = part_val + (depth-1)*r;
					for(int i=0 ; i< r; i++)
					{
						target[i] = matvali[i]*source[i];
					}
				}
				depth --;
				// Intermediate values (KrP) are calculated
				// and now written back to intval
				/*
				const TYPE* krp_res = part_val + (depth) * r;
				TYPE  * const intval_loc = t->intval[mode-1] + pos*r;
				for(int i=0 ; i< r; i++)
				{
					intval_loc[i] = krp_res[i];
					//cout<<krp_res[i]<<endl;
				}
				*/
				
				pos ++;
				
				for(idx_t f = t->ptr[mode-1][stack[depth]] ; f < t->ptr[mode-1][stack[depth]+1] ; f++)
				{

					idx_t fid = t->ind[mode][f];
					TYPE * const target = mats[mode]->val + ((mats[mode]) -> dim2) * fid; 
					const TYPE val = t->val[f];
					const TYPE* pp = part_val + (mode-1)*r;
					//cout<<stack[depth-1]<<" "<<stack[depth] <<" "<<fid<<" "<<depth<<endl;
					for(int i=0; i<r ; i++)
					{
						target[i] += val*pp[i];
					}
				}
				
				stack[depth] ++;
				while(depth >0 && t->ptr[depth-1][stack[depth-1]+1] == stack[depth])
				{
					depth --;
					stack[depth]++;
				}
				
			}
		}
	}
	
	return 0;
}


int saxpy_krp(csf * t, matrix** mats,int r)
{
	const idx_t nmodes = t->nmode;
	idx_t mode = nmodes-1;
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		#ifdef OMP
		int th = omp_get_thread_num();
		#else
		int th = 0;
		#endif
		idx_t start  = t->thread_start[th];
		idx_t stop = t->thread_start[th+1];
		TYPE* const part_val = new TYPE[nmodes * r];
		idx_t* stack = new idx_t[nmodes];
		idx_t depth = 0;
		for(idx_t s=start ; s<stop ; s++)
		{
			stack[0] = s;
			idx_t sid = t->ind[0][s];
			const TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * sid; 
			
			for(int i=0 ; i< r; i++)
			{
				part_val[i] = matval[i];
			}
			depth ++;
			idx_t pos=s;
			for(int i=0; i<nmodes-2 ; i++)
			{
				pos = t->ptr[i][pos];
			}
			for(int i=1; i<nmodes ; i++)
			{
				stack[i] = t->ptr[i-1][stack[i-1]];
				
			}
			
			
			while( depth > 0 && stack[0] == s)
			{
				for(;depth < mode; depth++)
				{
					idx_t sfid = t->ind[depth][stack[depth]];
					const TYPE* matvali = mats[depth]->val + ((mats[depth]) -> dim2) * sfid; 
					TYPE * const target = part_val + depth*r;
					const TYPE* source = part_val + (depth-1)*r;
					for(int i=0 ; i< r; i++)
					{
						target[i] = matvali[i]*source[i];
					}
				}
				depth --;
				// Intermediate values (KrP) are calculated
				// and now written back to intval
				const TYPE* krp_res = part_val + (depth) * r;
				TYPE  * const intval_loc = t->intval[mode-1] + pos*r;
				for(int i=0 ; i< r; i++)
				{
					intval_loc[i] = krp_res[i];
					//cout<<krp_res[i]<<endl;
				}
				
				
				pos ++;
				
				
				stack[depth] ++;
				while(depth >0 && t->ptr[depth-1][stack[depth-1]+1] == stack[depth])
				{
					depth --;
					stack[depth]++;
				}
				
			}
		}
	}
	return 0;
}

int saxpy_reduce(csf * t, matrix** mats, int r, Reverse_Indices* ri)
{
	int nmode = t->nmode;
	int mode = nmode - 1;
	
	memset(mats[mode]->val, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));
	LIKWID_MARKER_INIT;
	
	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		{
			LIKWID_MARKER_THREADINIT;	
		}
	}	

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
	//	if(profile == mode)
		{
			LIKWID_MARKER_START("Compute");
		}
	}	

	#ifdef OMP
	#pragma omp parallel
	#endif
	{
		#pragma omp for	schedule(guided,1024)
		for(idx_t row=0; row< t->mlen[mode] ; row++ )
		{
			TYPE * __restrict__ const matval = mats[mode]->val + ((mats[mode]) -> dim2) * row;
			for(idx_t nnz=ri->ptr[row] ; nnz < ri->ptr[row+1] ; nnz++)
			{
				const TYPE * __restrict__ const krp = t->intval[mode-1] + ri->krp[nnz]*r;
				const TYPE val = *(ri->vals[nnz]);
			//cout<<ri->krp[nnz]<<" "<<(ri->vals[nnz]) <<endl;
				#pragma omp simd
				for(int i=0 ; i< r ; i++)
				{
					matval[i] += krp[i]*val;
				}
			}
		}
		{
			LIKWID_MARKER_STOP("Compute");
		}
	}
	LIKWID_MARKER_CLOSE;
	return 0;
	
}


int main(int argc, char** argv)
{


	if(argc == 1)
	{
		printf("Usage is %s <matrix name> <number of ranks (optional)> \
<kernel (-1: SpTL , -2: SpTLso , 0-d!: permutation of modes) (optional)> \
<mode to profile (optional) >  \n", argv[0] );
		exit(1);
	}

	int nmode,i,r,mode;
	int debug = 1;
	csf* t = malloc_csf();
	coo* dt = NULL; 
	int profile = -1;
	int order_num = -1;
	if (argc > 3)
		order_num = atoi(argv[3]);
	if (argc > 4)
		profile = atoi(argv[4]);
	if(debug)
	{
		dt = malloc_coo();
		read_tensor(argv[1],t,dt,order_num);
	}
	else
	{
		read_tensor(argv[1],t,dt,order_num);
	}

	t->intval = NULL;

	
	print_csf(t,argv[1]);
	
	matrix** mats;
	nmode = t->nmode;


	count_fiber_leaf_root_fast(t);

	r = 32;
	if(argc > 2)
		r = atoi(argv[2]);

	mats = (matrix **) malloc(nmode*sizeof(matrix*));
	for(i=0 ; i<nmode ; i++)
	{
		mats[i] = create_matrix(t->mlen[i],r,1);
	}
	for(i=0 ; i<nmode ; i++)
	{
		random_matrix(*mats[i],i);
		if(VERBOSE  == VERBOSE_DEBUG)
			print_matrix(*mats[i]);
	}

	mttkrp_fused_init(t,r,true);

	double total=0;

//	printf("here\n")
	/*
	for(mode = 0 ; mode<nmode ; mode++)
	{
		mttkrp_atomic(t,mode,r,mats,profile);
		random_matrix(*mats[mode],i);
	}
	*/
	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	/*
	{
		auto start = std::chrono::high_resolution_clock::now();
		mttkrp_hardwired_first(t,0,r,mats,profile);
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		//total += diff.count();	
		printf("Hardwired time for mode %d %lf \n",t->modeid[0],diff.count());	
		//mttkrp_test(dt,0,r,mats);
	}
	*/

	bool correctness_error = false;

	double* times = new double[nmode];

	for(mode = nmode-1 ; mode<nmode ; mode++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		bool is_atomic = mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) );

		if(is_atomic)
			mttkrp_atomic(t,mode,r,mats,profile);
		else
			mttkrp_hardwired(t,mode,r,mats,profile);
		//printf("here %lf\n",mats[mode]->val[0] );
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		total += diff.count();
		if (is_atomic)
			printf("atomic     ");
		else
			printf("privatized ");

		printf("time for mode %d %lf \n",t->modeid[mode],diff.count());
		times[t->modeid[mode]] = diff.count();

		if(debug)
		{

			auto start2 = std::chrono::high_resolution_clock::now();
			int num_diff = mttkrp_test(dt,mode,r,mats);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end2-start2;
			//total += diff.count();
			if(num_diff)
			{
				correctness_error = true;
				//printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}
		}
		random_matrix(*mats[mode],i);

		for(i=0 ; i<nmode ; i++)
		{
			
			if(VERBOSE  == VERBOSE_DEBUG)
			{
				print_matrix(*mats[i]);
				
			}
			//random_matrix(*mats[i],i+1);
			//set_matrix(*mats[i],1);
		}


	}

	if (correctness_error)
	{
		printf("Results of MTTKRP is not correct\n");
		printf("incorrect,");
	}	

	if(order_num == -1)
		printf("SpTL1,");
	if(order_num == -2)
		printf("SpTL1so,");
	printf("%s", argv[1]);
	for(mode = 0 ; mode<nmode ; mode++)
		printf(",%lf", times[mode]);

	printf("\n");

	printf("Total MTTKRP time %lf \n",total);

	{
		auto start = std::chrono::high_resolution_clock::now();
		saxpy(t,mats,r);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		printf("saxpy %lf \n",diff.count());
		if(debug)
		{

			auto start2 = std::chrono::high_resolution_clock::now();
			int num_diff = mttkrp_test(dt,t->nmode-1,r,mats);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end2-start2;
			//total += diff.count();
			if(num_diff)
			{
				correctness_error = true;
				//printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}
		}
	}

	Reverse_Indices* ri;
	
	{
		auto start = std::chrono::high_resolution_clock::now();
		ri = find_reverse_indices(t);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		printf("time reverse index computation %lf \n",diff.count());
	}
		
	
	{
		auto start = std::chrono::high_resolution_clock::now();
		saxpy_krp(t,mats,r);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		printf("saxpy krp %lf \n",diff.count());
	}

	{
		auto start = std::chrono::high_resolution_clock::now();
		saxpy_reduce(t,mats,r,ri);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		printf("saxpy reduce %lf \n",diff.count());
		if(debug)
		{

			auto start2 = std::chrono::high_resolution_clock::now();
			int num_diff = mttkrp_test(dt,t->nmode-1,r,mats);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end2-start2;
			//total += diff.count();
			if(num_diff)
			{
				correctness_error = true;
				//printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}
			else
			{
				cout<<"saxpy reduce correctness check passed"<<endl;
			}
		}
	}		
	
	

	free_csf(t);
	if(debug)
		free_coo(dt);
	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}
