#ifndef READER_CC
#define READER_CC

#include "../inc/reader.h"

int num_mode;
int sort_shift;
int* sort_order;

int create_perm(int pern_num, int* order, int size)
{
	int it = size;
	int n = size;
	int i,j;
	int remainder;
	int avail = 1;
	int next = pern_num;
	
	while(it > 0)
	{
		it --;
		remainder = next % n;
		next = next / n;
		//printf("remainder %d\n",remainder);
		i = 0;
		while(0 <= remainder)
		{
			avail = 1;
			for(j = n; j < size; j++)
			{
				if(order[j] == i)
					avail = 0;
			}
			if(avail)				remainder --;
			i ++;
		}

		n --;
		order[it] = i-1;
	}
	/*
	for(i = 0; i< size; i++)
	{
		printf("%d ",order[i] );
	}
	printf("is the order\n");
	*/
	return 0;
}

int compare_nnz(const void *a, const void *b)
{
	idx_t** x = (idx_t**) a;
	idx_t** y = (idx_t**) b;
	int shift = sort_shift;
	int i, ii;
	int size = num_mode ;//*(x+1) - *x;
	//printf(" size %d\n", size);
	for(ii=0 ; ii<size ; ii++)
	{
		//i = (ii+shift) % size;
		i = sort_order[ii];
		if((*x)[i] < (*y)[i])
			return -1;
		else if ((*x)[i] > (*y)[i])
			return 1;
	}
	return 1;
}

int count_fiber(idx_t** pindex, int nnz, int nmode, int shift, int* fiber_count, int* sort_order)
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

int printt(idx_t** pindex, int nnz, int nmode, TYPE* vals)
{
	int i,j;
	for(i = 0; i< nnz ; i++)
	{
		for(j =0;j<nmode;j++)
		{
			printf("%d ",pindex[i][j]);
		}
		printf("\n");
	}
	return 0 ;
}

int readline(char* line, idx_t* idx, TYPE* val)
{
	char* start = line, *end;
	int m_id = 0;
	char convert[50];
	end = convert;
	int isfloat = 0;
	int valinit = 1;
	*val = 0;
	while(*(start) != '\0')
	{
		if(isdigit(*start))
		{
			*end = *start;
			end ++;
		}
		else if(*start == '.')
		{
			*end = *start;
			end ++;
			isfloat = 1;
		}
		else
		{
			*end= '\0';
			if (isfloat == 0 || m_id == MAX_MODE)
			{
				idx[m_id] = atoi(convert);

				//printf("%d ",idx[m_id]);
				m_id ++;

			}
			else
			{
				*val = atof(convert);
				//printf("%lf %ld ",*val,start-line);
				valinit = 0;
			}
			end = convert;
		}
		start ++;
	}




	if (valinit) // last value of the read was the value of the nnz
	{
		m_id --;
		*val = (TYPE) idx[m_id];
	}

	return m_id;
}

int print_stats(int* sort_order, int* fiber_count, const char* file, int nmode, int nnz)
{
	int i, ii;
	
	printf("%s Active number of fibers for modes ",file);
	for(ii = 0; ii < nmode; ii++)
	{
		printf("%d ",sort_order[ii]);	
	}
	


	printf("is ");

	for(ii = 0; ii < nmode ; ii++)
	{
		printf("%d ",fiber_count[ii]);
	}
	printf("\n");


	printf("%s NNZ per fiber for mode order ",file);
	for(ii = 0; ii < nmode; ii++)
	{
		printf("%d ",sort_order[ii]);	
	}
	
	printf("is ");

	for(ii = 0; ii < nmode; ii++)
	{
		printf("%lf ", ((double) nnz)/fiber_count[ii]);
	}
	printf("\n");

	return 0;
}

int reorder_stat(int nmode, int nnz, idx_t** pindex, const char* file)
{
	int i,ii;
	int nfac = 1;
	int * fiber_count;
	for(i = 1; i <= nmode; i++)
		nfac *= i;

	fiber_count = (int*) malloc(nmode*sizeof(int));
	for(i = 1; i <= nfac; i++)
	{

		create_perm(i , sort_order , nmode);
		qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);

		int num_fiber = count_fiber(pindex,nnz,nmode,ii,fiber_count,sort_order);
		
		print_stats(sort_order, fiber_count, file, nmode, nnz);
			
	}
	return 0;
}

int order_modes(int* mlen, int nmode, int* sort_order)
{
	// Simple bubble sort of the mlem array and the result written in sort_order
	int i, ii, ind, min, *len;

	int *sorted = (int* ) malloc(nmode*sizeof(int)) ;
	len = (int* ) malloc(nmode*sizeof(int)) ;

	ii = 0;
	for(i = 0 ; i < nmode ; i++)
	{
		sort_order[i] = -1;
		sorted[i] = -1;
	}

	while (ii < nmode)
	{
		min = -1;
		for (i = 0 ; i<nmode ; i++)
		{
			if (sorted[i] == -1 && ( min == -1 || min > mlen[i]))
			{
				min = mlen[i];
				ind = i;
			}
		}

		sorted[ind] = ii;
		sort_order[ii] = ind;
		len[ii] = mlen[ind];
		ii ++;
	}

	for(i = 0 ;i<nmode; i++)
	{
		mlen[i] = len[i];
	}

	rem(sorted);
	rem(len);
	return 0;
}


int coo2csr(idx_t** pindex, idx_t* index, TYPE* vals, int nnz, int nmode, int* fiber_count, csf* res,int* mlen)
{
	csf t;
	int i, j, jj , ii, ilen, plen, len, *ind, *dimlen;
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
	t.fiber_count = (int* ) malloc(nmode*sizeof(int));

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

	printf("Total space requirement of the CSF is %ld%s \n",total_space,space_sign);

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


	for( i=0 ; i<nmode ; i++)
	{
		t.inds[t.ptr[i][0]] = pindex[0][i];
		ind[i] = 1;
		dimlen[i] = 1;
	}

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
			if(diff > 0)
			{
				int loc = t.ptr[jj][ind[jj]];
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

	t.mlen = (int*) malloc(nmode*sizeof(int));
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

int read_tensor(const char* file, csf* res)
{
	FILE *fp;
	int *loc;
	char buf[300];
	int nflag, sflag;
	int dummy, pre_count=0, tmp_ne;
	int i,j,ii;
	idx_t* index = NULL;
	idx_t** pindex = NULL;
	int sizestep = 1000000;
	int size = sizestep;
	int nmode = 0;
	int nnz = 0;
	int* mlen;
	int * fiber_count;
	TYPE* vals;
	srand(time(NULL));
	//int compare_nnz(const void *a, const void *b);


	fp = fopen(file, "r");

	//for(int i = 0 ; i < 10 ; i++)
	while(fgets(buf, (MAX_MODE+1)*30, fp))
	{
		if(strstr(buf, "#") != NULL)
			continue;
		idx_t idx[MAX_MODE];
		TYPE val;
		//printf("%s",buf );
		int mode_len = readline(buf,idx,&val);
		//printf("\n%d %lf \n",mode_len,val);
		if(index == NULL || pindex == NULL)
		{
			nmode = mode_len;
			num_mode = nmode;
			index = (idx_t*) malloc(size*nmode*sizeof(idx_t));
			vals = (TYPE*) malloc(size*sizeof(TYPE));
			pindex =(idx_t**)  malloc((size+1)*sizeof(idx_t*));
			mlen = (int*) malloc(nmode*sizeof(int));
			for(i = 0; i<nmode ; i++)
				mlen[i] = 0;
			//pindex[0] = index;
			
		}
		else if(nnz == size)
		{
			size += sizestep;
			index = (idx_t*) realloc(index,size*nmode*sizeof(idx_t));
			pindex = (idx_t**) realloc(pindex,(size+1)*sizeof(idx_t*));
			vals = (TYPE*) realloc(vals,size*sizeof(TYPE));

			if(index == NULL || pindex == NULL || vals == NULL)
			{
				printf("Insufficient memory when reading");
				exit(1);
			}
		}

		for(i = 0; i<nmode ; i++)
		{
			index[nnz*nmode + i] = idx[i];
			if(idx[i] + 1  > mlen[i])
			{
				mlen[i] = idx[i]+1;
			}
		}
		vals[nnz] = val;

		nnz ++;
		//pindex[nnz] = index + nnz*nmode;
	} // Reading is complete


	fclose(fp);
	for(i = 0; i<size+1 ; i++)
		pindex[i] = index + nmode*i;

	printf("Tensor has dimensions ");
	for(i = 0; i<nmode-1 ; i++)
		printf("%dx",mlen[i]);
	printf("%d and %d nnz\n",mlen[i],nnz);

	sort_order = (int*) malloc(nmode*sizeof(int));
	fiber_count = (int*) malloc(nmode*sizeof(int));

	order_modes(mlen, nmode, sort_order);
	/*
	sort_order[0] = 0;
	sort_order[1] = 1;
	sort_order[2] = 2;
	*/
	qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);
	count_fiber(pindex,nnz,nmode,ii,fiber_count,sort_order);

	if(VERBOSE > VERBOSE_HIGH)
		print_stats(sort_order, fiber_count, file, nmode, nnz);


	//reorder_stat(nmode, nnz,  pindex, file );

	csf *t = res;
	coo2csr(pindex, index, vals, nnz,nmode,fiber_count,t,mlen);
	if(VERBOSE == VERBOSE_DEBUG)
	{
		printt(pindex , nnz , nmode, vals);

		printf("nmode %d \n ",t->nmode);
	}



	rem(index);
	rem(pindex);
	rem(vals);
	rem(sort_order);
	rem(mlen)
	rem(fiber_count);

	return 0;
}



#endif