#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

int main(int argc, char** argv)
{
	int nmode,i,r,mode;
	csf* t = (csf *) malloc(sizeof(csf));
	read_tensor(argv[1],t);
	matrix** mats;
	nmode = t->nmode;

	r = 32;
	if(argc > 2)
		r = atoi(argv[2]);

	mats = (matrix **) malloc(nmode*sizeof(matrix*));
	for(i=0 ; i<nmode ; i++)
	{
		mats[i] = create_matrix(t->mlen[i],r,i+1);
	}
	for(i=0 ; i<nmode ; i++)
	{
		if(VERBOSE  == VERBOSE_DEBUG)
			print_matrix(*mats[i]);
	}

	mttkrp_fused_init(t,r);

	for(mode = 0 ; mode<nmode ; mode++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		mttkrp_atomic(t,mode,r,mats);
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;

		printf("time for mode %d %lf \n",t->modeid[mode],diff.count());
		for(i=0 ; i<nmode ; i++)
		{
			if(VERBOSE  == VERBOSE_DEBUG)
				print_matrix(*mats[i]);
		}
	}


	
	

	//free(t->intval[1]);
	free_csf(t);

	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}