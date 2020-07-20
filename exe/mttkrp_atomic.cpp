#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

int main(int argc, char** argv)
{
	int nmode,i,r,mode;
	int debug = 1;
	csf* t = (csf *) malloc(sizeof(csf));
	coo* dt = NULL; 
	if(debug)
	{
		dt = (coo *) malloc(sizeof(coo));
		read_tensor(argv[1],t,dt);
	}
	else
	{
		read_tensor(argv[1],t);
	}

	t->intval = NULL;
	int profile = -1;

	if (argc > 3)
		profile = atoi(argv[3]);

	print_csf(t);
	
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

//	printf("here\n")

	for(mode = 0 ; mode<nmode ; mode++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		mttkrp_atomic(t,mode,r,mats,profile);
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

	printf("Total MTTKRP time %lf \n",total);

	
	

	//free(t->intval[1]);
	free_csf(t);
	free_coo(dt);
	rem(t);
	rem(dt);
	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}