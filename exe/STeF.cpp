#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/mttkrp_combined.h"
#include "../inc/mttkrp_hardwired.h"
#include <time.h>
#include <string>
#include <cmath>  // For the pow function
// #define MAX(X, Y) X>Y ? X : Y


long estimate_data_movement(long cache_size, csf* t, int r, bool* intv)
{
	int nmode = t->nmode;	
	long no_memo_data = 0;
	for(int i=0 ; i < nmode; i++)
	{
		no_memo_data += (long) 2*t->fiber_count[i];
		long mat_size = (long) t->mlen[i] * r;
		no_memo_data += mat_size > cache_size ? (long) t->fiber_count[i] * r : mat_size;
	}

	bool some_memoized = false;
	for(int i = 0 ; i < nmode-2; i++)
		if (intv[i])
			some_memoized = true;
	
	if (!some_memoized)	
		return no_memo_data*nmode;

	long res = 0;	
	for(int i=0 ; i < nmode; i++)
	{
		int memo_idx = nmode;
		for(int j=MAX(0, i-1); j<nmode-2; j++)
			if(intv[j])
			{
				memo_idx = j+1;
				break;
			}
		
		if (memo_idx == nmode)
		{
			res += no_memo_data;
		}
		else
		{
			long partial = 0;
			for(int j=0; j<memo_idx; j++)
			{
				partial += (long) 2*t->fiber_count[j];
				long mat_size = (long) t->mlen[j] * r;
				partial += mat_size > cache_size ? (long) t->fiber_count[j] * r : mat_size;
				partial += mat_size;
			}
			res += partial;
		}		
	}
	return res;
}

void covert_int_to_bool_arr(bool* intv, int case_idx, int len)
{
	int rem = case_idx;
	for(int i = 0; i < len; i++)
	{
		intv[i] = (rem % 2 == 0);
		rem /= 2;
	}
}

bool* get_model(long cache_size, csf* t, int r)
{
	int nmode = t->nmode;
	bool* intv = new bool[nmode-2];
	long best_score = 0;
	int best_case_idx = 0;
	for (int case_idx = 0 ; case_idx < (1 << (nmode-2)) ; case_idx ++)
    {
		covert_int_to_bool_arr(intv, case_idx, nmode-2);
		long case_score = estimate_data_movement(cache_size, t, r, intv);
		if (case_idx == 0 || case_score < best_score)
		{
			best_score = case_score;
			best_case_idx = case_idx;
		}	
		printf("Configuration %d estimated data movement %ld\n",case_idx, case_score);
	}
	printf("Configuration %d chosen by model. Estimated data movement %ld\n",best_case_idx, best_score);
	covert_int_to_bool_arr(intv, best_case_idx, nmode-2);
	return intv;
}


int main(int argc, char** argv)
{
	int nmode,i,r;
	int debug = 1;
	csf* t = malloc_csf();
	coo* dt = NULL; 
	int profile = -1;
	int order_num = -1;

	if (argc < 2)
	{
		printf("Usage is %s <tensor> <number of ranks>\n",argv[0]);
		return 0;
	}
	if(debug)
	{
		dt = malloc_coo();
		read_tensor(argv[1],t,dt,order_num);
	}
	else
	{
		read_tensor(argv[1],t,dt,order_num);
	}
	
	r = 32;
	if(argc > 2)
		r = atoi(argv[2]);

	long cache_size = 35*1024*1024;
	if(argc > 3)
		cache_size = atoi(argv[3]);

	t->intval = NULL;

	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	print_csf(t,argv[1]);
	
	matrix** mats;
	nmode = t->nmode;

	mats = (matrix **) malloc(nmode*sizeof(matrix*));
	for(i=0 ; i<nmode ; i++)
	{
		if(i == 0)
			mats[i] = create_matrix(t->mlen[i],r,1);
		else
			mats[i] = create_matrix(t->mlen[i],r,1);
	}
	for(i=0 ; i<nmode ; i++)
	{
		random_matrix(*mats[i],i);
		if(VERBOSE  == VERBOSE_DEBUG)
			print_matrix(*mats[i]);
	}

	//mttkrp_fused_init(t,r,true);

	double total=0;
	b_thread_start(t);
	mttkrp_fused_init_ms(t,r,true,NULL);
	
	if (nmode > 5)
	{
		printf("Only supported tensor with at most 5 dimensions.\n");
		return 0;
	}

    

	// Find the model
    // for (int repeat = 0 ; repeat < num_cases ; repeat ++)
    // {
    //     bool* intv = new bool[nmode-2];
    //     {
    //         int rem = repeat;
    //         for(int i = 0; i < nmode-2 ; i++)
    //         {
    //             intv[i] = (rem % 2 == 0);
    //             rem /= 2;
    //         }
    //     }
	{
		bool* intv = get_model(cache_size, t, r);
        char* save_str = new char[nmode-1];
        for(int i = 0; i < nmode-2 ; i++)
            if (intv[i])
                save_str[i] = 's';
            else
                save_str[i] = 'n';

        save_str[nmode-2] = '\0';
        total = 0;
        for(int mode = 0 ; mode<nmode ; mode++)
        {
            // clock_t cstart, cend;
            bool is_atomic = (mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) ))	;
            auto start = std::chrono::high_resolution_clock::now();
            // cstart = clock();
            mttkrp_combined_lb(t,r,mats,profile,mode,!is_atomic,intv);            
            // cend = clock();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end-start;
            total += diff.count();
            printf("IS is %s, %s for mode %d %lf \n",save_str,(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
            // double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
            // printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
            if(debug)
            {
                //  for(i=0 ; i<nmode ; i++)  if (i != mode)               random_matrix(*mats[i],i);
                auto start2 = std::chrono::high_resolution_clock::now();
                mttkrp_test(dt,mode,r,mats);
                auto end2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff = end2-start2;
                //total += diff.count();
                //printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
            }
            for(i=0 ; i<nmode ; i++)
            {
                
                if(VERBOSE  == VERBOSE_DEBUG)
                {
                    print_matrix(*mats[i]);
                    
                }
                //random_matrix(*mats[i],i+1);
                //set_matrix(*mats[i],1);
            }
            random_matrix(*mats[mode],mode);
        }
        printf("Total Intermediate Save %s combined time template MTTKRP time %lf\n",save_str,total);
		// 0 intval caches
		
		for(int i=1;i<nmode-1;i++)
		{
			if(t->intval != NULL && intv[i-1])
				memset(t->intval[i],0,sizeof(TYPE)*(t->fiber_count[i] + t->num_th)*r);
		}
    }

	//free(t->intval[1]);
	free_csf(t);
	free_coo(dt);
	//rem(t);
	//rem(dt);
	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}
