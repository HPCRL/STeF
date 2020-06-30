#ifndef READER_CC
#define READER_CC

#include "../inc/reader.h"

int num_mode;
int sort_shift;

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
		i = (ii+shift) % size;
		if((*x)[i] < (*y)[i])
			return -1;
		else if ((*x)[i] > (*y)[i])
			return 1;
	}
	return 1;
}

int count_fiber(idx_t** pindex, int nnz, int nmode, int shift)
{
	int num_fiber = nnz;
	int i,j,jj,diff;
	
	for(i = 1; i<nnz; i++)
	{
		diff = 0;
		for(jj=0;jj<nmode-1;jj++)
		{
			j = (jj+shift) % nmode;
			if(pindex[i][j] != pindex[i-1][j])
				diff ++;
		}
		if (diff == 0)
		{
			num_fiber --;
		}
	}
	return num_fiber;
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
			pindex[0] = index;
			
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
		}

		nnz ++;
		pindex[nnz] = index + nnz*nmode;
	} // Reading is complete

	for(ii = 0; ii < nmode; ii++)
	{
		//for(i = 0; i<nnz ; i++)
		//	printf("%ld ", (pindex[i] - pindex[0])/3 + 1);
		//printf("\n");
		sort_shift = ii;
		qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);
		//for(i = 0; i<nnz ; i++)
		//	printf("%ld ", (pindex[i] - pindex[0])/3 + 1);
		//printf("\n");
		int num_fiber = count_fiber(pindex,nnz,nmode,ii);
		double nnz_per_fiber = ((double) nnz)/num_fiber;
		printf("Active number of fibers for mode %d is %d\n",(nmode - 1 + ii)%nmode +1, num_fiber);
		printf("NNZ per fiber for mode %d is %lf\n",(nmode - 1 + ii)%nmode +1, nnz_per_fiber);
	}
	//if(strstr(buf, "symmetric") != NULL || strstr(buf, "Hermitian") != NULL) sflag = 1;

	return 0;
}



#endif