#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

int main(int argc, char** argv)
{
	int nmode,i,r,mode;
	int debug = 0;
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

	
	print_csf(t);
	
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

	mttkrp_fused_init(t,r);

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
	for(mode = 0 ; mode<nmode ; mode++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		if(mode == 0 || t->fiber_count[mode] / t->mlen[mode] > num_th * ATOMIC_THRESH )
			mttkrp_hardwired(t,mode,r,mats,profile);
		else
			mttkrp_atomic(t,mode,r,mats,profile);
		//printf("here %lf\n",mats[mode]->val[0] );
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		total += diff.count();
		printf("time for mode %d %lf \n",t->modeid[mode],diff.count());

		if(debug)
		{
			auto start2 = std::chrono::high_resolution_clock::now();
			mttkrp_test(dt,mode,r,mats);
			auto end2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end2-start2;
			//total += diff.count();

			printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
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


	/*
	{

		auto start = std::chrono::high_resolution_clock::now();
		mttkrp_hardwired_last(t,nmode-1,r,mats,profile);
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		
		printf("Hardwired time for mode %d %lf \n",t->modeid[nmode-1],diff.count());	
		mttkrp_test(dt,nmode-1,r,mats);
	}
	*/
	printf("Total MTTKRP time %lf \n",total);

	
	

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
