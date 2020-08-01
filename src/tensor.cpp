#ifndef TENSOR_CPP
#define TENSOR_CPP


#include "../inc/tensor.h"



csf* malloc_csf()
{
	csf* res = (csf*) malloc(sizeof(csf));
	res->ptr = NULL;
	res->ptrs = NULL;
	res->ind = NULL;
	res->inds = NULL;
	res->val = NULL;
	res->fiber_count = NULL;
	res->mlen = NULL;
	res->modeid = NULL;
	res->intval = NULL;
	return res;
}

coo* malloc_coo()
{
	coo* res = (coo*) malloc(sizeof(coo));
	res -> ind = NULL;
	res -> val = NULL;
	return res;
}


int free_csf(csf* t)
{
	int i;
	rem(t->ptr);
	rem(t->ptrs);
	rem(t->ind);
	rem(t->inds);
	rem(t->val);
	rem(t->fiber_count);
	rem(t->mlen);
	rem(t->modeid);
	//printf("%d\n", t->nmode );
	if(t->intval != NULL)
		for (i=0; i< (t->nmode) ; i++)
			rem(t->intval[i]);
	if(t->private_mats != NULL)
		for (i=0; i< (t->num_th) ; i++)
			free_matrix(t->private_mats[i]);
	rem(t->private_mats);
	rem(t->intval);
	rem(t);
	return 0;
}

int free_coo(coo* t)
{
	rem(t->ind);
	rem(t->val);
	rem(t);
	return 0;
}


int print_csf(csf* t)
{
	for(int i=0; i<t->nmode ;i++)
	{
		printf("Mode %d fiber count %d\n", i, t->fiber_count[i]);
		printf("Mode %d mode length %d\n", i, t->mlen[i]);
		//printf("Mode %d fiber_count %d\n", fiber_count[i]);
	}
	return 0;
}

int csf_space(csf* t)
{

	return 0;
}

int find_inds(idx_t* inds ,csf* t,idx_t it)
{
	int nmode = t->nmode;
	//printf("here ss\n");
	if (it == 0)
	{
		for (int i=0; i<nmode ; i++)
		{
			inds[i] = 0;
		}

		return 0;
	}
	else
	{
		
		inds[nmode-1] = it;
		idx_t last_pos = it;
		for(int i=nmode-2 ; i>=0 ; i--)
		{
			// do binary seach and find the index
			// idx_t id = -1;
			idx_t start = 0;
			idx_t end = t->fiber_count[i];
			while(end > start+1)
			{
				idx_t pivot = (end + start)/2;
				if(t->ptr[i][pivot] > last_pos)
				{
					end = pivot;
				}
				else if (t->ptr[i][pivot] <= last_pos)
				{
					start = pivot;
				}
				//printf("%d %d %d %d\n",start, end , i , last_pos);
			}
			inds[i] = start;
			last_pos = start;
		}
		/*
		for(int i=0; i<nmode ; i++)
		{
			printf("%d", inds[i] );
		}
		printf("\n");
		*/
		return 0;
	}

	
}



int coo2csf(idx_t** pindex, idx_t* index, TYPE* vals, idx_t nnz, int nmode, idx_t* fiber_count, csf* res,int* mlen, int* sort_order)
{
	csf t;
	int i, j, jj , ii, ilen, plen,  *ind, *dimlen;
	long long total_space;
	char* space_sign;

	t = *res;
	plen = 0;
	for(i = 0; i < nmode -1  ; i++)
		plen += fiber_count[i] + 1 ;

	ilen = plen +  nnz + 1;
	plen ++;
	t.inds = (idx_t* ) malloc(ilen*sizeof(idx_t));
	t.ptr = (idx_t** ) malloc((nmode+1)*sizeof(idx_t*));
	t.ind = (idx_t** ) malloc((nmode+1)*sizeof(idx_t*));
	t.ptrs = (idx_t* ) malloc((plen)*sizeof(idx_t));
	t.val = (TYPE* ) malloc(nnz*sizeof(TYPE));
	t.modeid = (int* ) malloc(nmode*sizeof(int));
	dimlen = (int* ) malloc(nmode*sizeof(int));
	ind = (int* ) malloc(nmode*sizeof(int));
	t.fiber_count = (idx_t* ) malloc(nmode*sizeof(idx_t));

	ii = 0;
	

	total_space = ilen*sizeof(idx_t);
	total_space += 2*(nmode+1)*sizeof(idx_t*);
	total_space += plen*sizeof(idx_t);
	total_space += nnz*sizeof(idx_t);
	total_space += 3*nmode*sizeof(idx_t);

	if(total_space >= 1073741824)
	{
		// GB
		space_sign = "GB";
		total_space /= 1073741824;
	}
	else if(total_space >= 1048576)
	{
		// MB
		space_sign = "MB";
		total_space /= 1048576;
	}
	else if(total_space >= 1024)
	{
		// KB
		space_sign = "KB";
		total_space /= 1024;
	}
	else
	{
		// B
		space_sign = "B";
	}

	printf("Total space requirement of the CSF is %llu%s \n",total_space,space_sign);

	for( i=0 ; i< plen ; i++)	
	{
		t.ptrs[i] = 0;
	}

	for( i=0 ; i<nmode-1 ; i++)
	{
		t.ptr[i] = t.ptrs + ii;
		t.ind[i] = t.inds + ii;
		ii += fiber_count[i] + 1;
		
	}
	t.ptr[nmode - 1]  = t.ptr[nmode] = t.ptrs + ii ;
	t.ind[nmode - 1] = t.inds + ii ;
	ii += nnz + 1;
	t.ind[nmode] = t.inds + ii ;
	/*
	printf("test %d\n",t.ind[0][0] );
	printf("%d %d %d %d \n", pindex[0][0], pindex[0][1], pindex[0][2], pindex[0][3]);
	printf("%d %d %d %d \n", sort_order[0], sort_order[1], sort_order[2], sort_order[3]);

	printf("%d %d %d %d \n", t.ptr[0][0], t.ptr[1][0], t.ptr[2][0], t.ptr[3][0]);
	*/
	for( i=0 ; i<nmode ; i++)
	{
		j = sort_order[i];
		t.ind[i][0] = pindex[0][j];
		ind[i] = 1;
		dimlen[i] = 1;
	//	printf("test %d\n",t.ind[0][0] );
	}



	//printf("test %d\n",t.ind[0][0] );

	for( i=1 ; i<nnz ; i++)
	{
		int diff = 0;
		for(jj = 0 ; jj < nmode  ; jj++)
		{
			j = sort_order[jj];
			if(pindex[i][j] != pindex[i-1][j] && diff == 0)
			{
				diff ++;
				if(jj > 0)
				{
					dimlen[jj-1] ++; 
				}
			}
			if(diff > 0 )
			{
			//	int loc = t.ptr[jj][ind[jj]];
				if(jj < nmode-1)
				{
					t.ptr[jj][ind[jj]] = t.ptr[jj][ind[jj]-1] + dimlen[jj];
					dimlen[jj] = 1;
				}
				t.ind[jj][ind[jj]] = pindex[i][j];
				ind[jj] ++;
			}
		}
	}
	for(i = 0 ; i<nmode - 1 ; i++)
	{
		t.ptr[i][ind[i]] = t.ptr[i][ind[i]-1] + dimlen[i];
	}

	for(i = 0 ; i<nmode  ; i++)
	{
		t.ind[sort_order[i]][ind[sort_order[i]]] = -1;
	}

	for(i = 0 ; i<nmode  ; i++)
	{
		t.fiber_count[i] = ind[i];
	}



	for(i = 0 ; i<nnz ; i++)
	{
		t.val[i] = *( vals  + (pindex[i] - index)/nmode);
	}

	if(VERBOSE == VERBOSE_DEBUG)
	{
		for(i = 0 ; i<nmode ; i++)
			printf("ind %d\n",ind[i]);
	
		printf("\n");
		printf("\n");
	
		for(i = 0 ; i<nmode+1 ; i++)
		{
			printf("%ld\n", t.ptr[i] - t.ptr[0]);
		}
	
		printf("\n");
	
		for(i = 0 ; i<nmode+1 ; i++)
		{
			printf("%ld\n", t.ind[i] - t.ind[0]);
		}
	
		printf("\n");
	
		for(i = 0 ; i<ilen ; i++)
			printf("%d %d %d %lf\n", i, t.ptrs[i], t.inds[i], t.val[i]);
		printf("\n");
	}
	t.nmode = nmode;

	t.mlen = (idx_t*) malloc(nmode*sizeof(idx_t));
	for(i = 0; i<nmode ; i++)
	{
		t.mlen[i] = mlen[i];
	}
	for(i = 0; i<nmode ; i++)
	{
		t.modeid[i] = sort_order[i];
	}


	*res = t;

	rem(dimlen);
	rem(ind);
	return 0;
}

int* tensor_sort_order;
int tensor_num_mode;


// Helper function for COO to CSF function

int tensor_compare_nnz(const void *a, const void *b)
{
	idx_t** x = (idx_t**) a;
	idx_t** y = (idx_t**) b;
	//int shift = sort_shift;
	int i, ii;
	int size = tensor_num_mode ;//*(x+1) - *x;
	//printf(" size %d\n", size);
	for(ii=0 ; ii<size ; ii++)
	{
		//i = (ii+shift) % size;
		i = tensor_sort_order[ii];
		if((*x)[i] < (*y)[i])
			return -1;
		else if ((*x)[i] > (*y)[i])
			return 1;
	}
	return 1;
}


int coo2csf(coo* dt, csf* t, int* sort_order)
{
	t = malloc_csf();
	idx_t nnz = dt->nnz;
	int nmode = dt->nmode;


	for(int i =0 ; i<nmode ; i ++)
		tensor_sort_order[i] = sort_order[i];

	idx_t** pindex = &(dt->ind);
	idx_t* fiber_count = (idx_t*) malloc(sizeof(idx_t)*nmode);
	int* mlen = (int*) malloc(sizeof(int)*nmode);

	qsort(pindex,nnz,sizeof(idx_t*),tensor_compare_nnz);
	count_fiber(pindex,nnz,nmode,-1,fiber_count,sort_order);

	for(int i =0 ; i<nmode ; i ++)
		mlen[i] = 1;

	for (int i=0 ;i<nnz ; i++)
	{
		for(int j =0 ; j<nmode ; j ++)
		{
			if(dt->ind[i*nmode + j] + 1 > mlen[j] )
			{
				mlen[j] = dt->ind[i*nmode + j] + 1;
			}		
		}
	}

	coo2csf(pindex,dt->ind,dt->val,nnz,nmode,fiber_count,t,mlen,sort_order);

	return 0;
}


int count_fiber(idx_t** pindex, idx_t nnz, int nmode, int shift, idx_t* fiber_count, int* sort_order)
{
	int num_fiber = nnz;
	int i,j,jj,diff;
	
	for(i = 0; i<nmode ; i++)
	{
		fiber_count[i] = 0;
	}

	for(i = 1; i<nnz; i++)
	{
		diff = 0;
		for(jj=0;jj<nmode-1;jj++)
		{
			//j = (jj+shift) % nmode;
			j = sort_order[jj];
			if(pindex[i][j] != pindex[i-1][j])
			{
				if(diff == 0)
				{
					fiber_count[jj] ++;
				}
				diff ++;
			}
		}
		if (diff == 0)
		{
			num_fiber --;
		}
	}
	fiber_count[0] ++;	
	/*
	printf("%d ",fiber_count[0] );
	*/
	for(i = 1; i < nmode - 1; i++)
	{
		fiber_count[i] += fiber_count[i-1];
	//	printf("%d ",fiber_count[i] );
	}
	fiber_count[nmode-1] = nnz;
	//printf("\n" );
	
	return num_fiber;
}


int count_fiber(coo* dt, int* sort_order, int hmode)
{
	int size = dt->nmode;
	idx_t nnz = dt->nnz;
	int nmode = dt->nmode;
	uint32_t hashkey = 0;
	uint32_t* tohash = new uint32_t[hmode];

	uint32_t hashres = hashword(tohash,size,hashkey);

	std::unordered_set<uint32_t>* sets = new std::unordered_set<uint32_t>[hmode];


	printf("hash res test is  %d\n", hashres);
	for(int i=0; i<nnz ; i++)
	{
		for(int m=0 ; m < hmode ; m++)
		{
			tohash[m] = dt->ind[i*nmode + sort_order[m]];
			hashres = hashword(tohash,m+1,hashkey);
			sets[m].insert(hashres);
		}
	}

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