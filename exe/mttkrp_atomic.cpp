#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

int main(int argc, char** argv)
{
	int nmode,i,r;
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

	mttkrp_atomic(t,nmode-1,r,mats);
	for(i=0 ; i<nmode ; i++)
	{
		if(VERBOSE  == VERBOSE_DEBUG)
			print_matrix(*mats[i]);
	}


	return 0;
}