#ifndef SAMPLE_CPP
#define SAMPLE_CPP

#include "../inc/sample.h"

int estimate_fiber(coo* dt, int* sort_order, int* fiber_count, double sample_rate)
{
	int size = dt->nmode;
	idx_t nnz = dt->nnz;
	int nmode = dt->nmode;
	uint32_t hashkey = 0;
	int hmode = nmode;
	uint32_t* tohash = new uint32_t[hmode];

	uint32_t hashres = hashword(tohash,size,hashkey);

	srand (time(NULL));

	std::unordered_set<uint32_t>* sets = new std::unordered_set<uint32_t>[hmode];

	std::unordered_set<uint32_t> sample;

	idx_t range = nnz * sample_rate;

	printf("hash res test is  %d\n", hashres);
	for(int ii=0; ii<range ; ii++)
	{
		uint32_t i = (uint32_t) rand() % nnz;
		while ( sample.find(i) != sample.end() )
		{
			// If the key is in set get a new random value
			i = (uint32_t) rand() % nnz;
		}
		for(int m=0 ; m < hmode ; m++)
		{
			tohash[m] = dt->ind[i*nmode + sort_order[m]];
			hashres = hashword(tohash,m+1,hashkey);
			sets[m].insert(hashres);
		}
	}

	printf("Estimate fiber results: \n");
	printf("For mode order of ");
	for(int i=0; i<hmode ; i++)
	{
		printf(" -> %d ", sort_order[i] );
	}
	printf(" number of fibers with hashmap is \n");

	for(int i=0; i<hmode ; i++)
	{
		printf("Mode %d, hash val, %ld \n", sort_order[i], sets[i].size() );
	}


	

	// someone sets size a positive value 
	//std::unordered_set<verylong*, MyHash, MyEqual> set(bucket_count, MyHash(size));
	return 0;
}

#endif