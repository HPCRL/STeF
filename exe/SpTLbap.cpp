#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/mttkrp_combined.h"
#include "../inc/mttkrp_hardwired.h"
#include <time.h>
#include <string>

int main(int argc, char** argv)
{
	int nmode,i,r;
	int debug = 1;
	csf* t = malloc_csf();
	coo* dt = NULL; 
	int profile = -1;
	int order_num = -1;
	if (argc > 3)
		order_num = atoi(argv[3]);
	if (argc > 4)
		profile = atoi(argv[4]);
	if (argc > 5)
		debug = atoi(argv[5]);
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

	int num_th = 1;
	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	print_csf(t,argv[1]);
	
	matrix** mats;
	nmode = t->nmode;

	r = 32;
	if(argc > 2)
		r = atoi(argv[2]);

	mats = (matrix **) malloc(nmode*sizeof(matrix*));
	for(i=0 ; i<nmode ; i++)
	{
		if(i == 0)
			mats[i] = create_matrix(t->mlen[i]+num_th,r,1);
		else
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

	

	b_thread_start(t);
	

	if (nmode > 5)
		return 0;

    int num_cases = 1;
    for(int i = 0 ; i< nmode-2; i++)
        num_cases *= 2;
    for (int repeat = 0 ; repeat < num_cases ; repeat ++)
    {
        bool* intv = new bool[nmode-2];
        {
            int rem = repeat;
            for(int i = 0; i < nmode-2 ; i++)
            {
                intv[i] = (rem % 2 == 0);
                rem /= 2;
            }
        }

        char* save_str = new char[nmode-1];
        for(int i = 0; i < nmode-2 ; i++)
            if (intv[i])
                save_str[i] = 's';
            else
                save_str[i] = 'n';

        save_str[nmode-2] = '\0';
        total = 0;
		for(int force_atomic = 0; force_atomic < 2 ; force_atomic++)
		{
			for(int mode = 0 ; mode<nmode ; mode++)
			{
				clock_t cstart, cend;
				bool is_atomic = (mode > 0 && (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) ) || (force_atomic == 1);
				auto start = std::chrono::high_resolution_clock::now();
				cstart = clock();
				mttkrp_combined_lb(t,r,mats,profile,mode,!is_atomic,intv);            
				cend = clock();
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end-start;
				total += diff.count();
				printf("IS is %s, %s for mode %d %lf \n",save_str,(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
				double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
				// printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
				if(debug == 1)
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
				if(intv[i-1])
					memset(t->intval[i],0,sizeof(TYPE)*(t->fiber_count[i] + t->num_th)*r);
			}
		}
    }
	

	
	

	//free(t->intval[1]);
	free_csf(t);
	if(debug)
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
