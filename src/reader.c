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
	//printf("\n" );
	
	return num_fiber;
}

int printt(idx_t** pindex, int nnz, int nmode)
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
			}
			end = convert;
		}
		start ++;
	}




	if (*val == 0) // last value of the read was the value of the nnz
	{
		m_id --;
		*val = (TYPE) idx[m_id];
	}

	return m_id;
}

int read_tensor(const char* file)
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
		}

		for(i = 0; i<nmode ; i++)
		{
			index[nnz*nmode + i] = idx[i];
			if(idx[i] + 1  > mlen[i])
			{
				mlen[i] = idx[i]+1;
			}
		}

		nnz ++;
		//pindex[nnz] = index + nnz*nmode;
	} // Reading is complete

	for(i = 0; i<size+1 ; i++)
		pindex[i] = index + nmode*i;

	printf("Tensor has dimensions ");
	for(i = 0; i<nmode-1 ; i++)
		printf("%dx",mlen[i]);
	printf("%d and %d nnz\n",mlen[i],nnz);

	int nfac = 1;
	for(i = 1; i <= nmode; i++)
		nfac *= i;

	sort_order = (int*) malloc(nmode*sizeof(int));
	fiber_count = (int*) malloc(nmode*sizeof(int));
	for(i = 1; i <= nfac; i++)
	{
		//for(i = 0; i<nnz ; i++)
		//	printf("%ld ", (pindex[i] - pindex[0])/3 + 1);
		//printf("\n");
		//sort_shift = ii;
		create_perm(i , sort_order , nmode);
		qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);

		//printt(pindex,nnz,nmode);
		//for(i = 0; i<nnz ; i++)
		//	printf("%ld ", (pindex[i] - pindex[0])/3 + 1);
		//printf("\n");
		int num_fiber = count_fiber(pindex,nnz,nmode,ii,fiber_count,sort_order);
		double nnz_per_fiber = ((double) nnz)/num_fiber;
		printf("%s Active number of fibers for modes ",file);
		for(ii = 0; ii < nmode; ii++)
		{
			printf("%d ",sort_order[ii]);	
		}
		
		printf("is ");

		for(ii = 0; ii < nmode - 1; ii++)
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

		for(ii = 0; ii < nmode - 1; ii++)
		{
			printf("%lf ", ((double) nnz)/fiber_count[ii]);
		}
		printf("\n");
			//%d is %lf\n",file,(nmode - 1 + ii)%nmode , nnz_per_fiber);
	}
	//if(strstr(buf, "symmetric") != NULL || strstr(buf, "Hermitian") != NULL) sflag = 1;

	//nmode = 5;
	
	
	return 0;
}



#endif