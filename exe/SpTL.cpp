#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

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

	for(mode = 0 ; mode<nmode ; mode++)
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
	printf("SpTL order=%d Total MTTKRP time %lf \n",order_num,total);

	
	

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
