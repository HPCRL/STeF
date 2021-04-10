#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/mttkrp_hardwired.h"
#include <time.h>

int main(int argc, char** argv)
{
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
	{

		auto start = std::chrono::high_resolution_clock::now();
		mttkrp_hardwired_first_not_fused(t,0,r,mats,profile);
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		
		printf("Hardwired time with no fusion for mode %d %lf \n",t->modeid[0],diff.count());	
		mttkrp_test(dt,0,r,mats);
	}

	/*
	for(mode = 0 ; mode<nmode ; mode++)
	{
		mttkrp_hardwired(t,mode,r,mats,profile);
		random_matrix(*mats[mode],i);
	}
	*/

	for(mode = nmode-1 ; mode<nmode ; mode++)
	{
		clock_t cstart, cend;
		auto start = std::chrono::high_resolution_clock::now();
		cstart = clock();
		mttkrp_hardwired(t,mode,r,mats,profile);
		cend = clock();
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		total += diff.count();
		printf("Hardwired time for mode %d %lf \n",t->modeid[mode],diff.count());
		double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
		printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
		if(debug)
		{
			auto start2 = std::chrono::high_resolution_clock::now();
			//mttkrp_test(dt,mode,r,mats);
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


	printf("Total Hardwired MTTKRP time %lf\n",total);
	

	
	

	free_csf(t);
	free_coo(dt);
	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}
