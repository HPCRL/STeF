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
	//int shift = sort_shift;
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
	while(*start != '\0')
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

	if(end != convert)
	{
		*end = '\0';
		*val = atof(convert);
		valinit = 0;
		end = convert;
	}

	if (valinit) // last value of the read was the value of the nnz
	{
		m_id --;
		*val = (TYPE) idx[m_id];
	}

	return m_id;
}

int print_stats(int* sort_order, idx_t* fiber_count, const char* file, int nmode, int nnz)
{
	int  ii;
	
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
	int i;
	int nfac = 1;
	idx_t * fiber_count;
	for(i = 1; i <= nmode; i++)
		nfac *= i;

	fiber_count = (idx_t*) malloc(nmode*sizeof(idx_t));
	for(i = 1; i <= nfac; i++)
	{

		create_perm(i , sort_order , nmode);
		qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);

		//int num_fiber = count_fiber(pindex,nnz,nmode,ii,fiber_count,sort_order);
		
		print_stats(sort_order, fiber_count, file, nmode, nnz);
			
	}
	return 0;
}


int order_modes(int* mlen, int nmode, int* sort_order)
{
	// Simple bubble sort of the mlem array and the result written in sort_order
	int i, ii, ind=-1, min, *len;

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


int sort_coo(idx_t** pindex,idx_t* index,TYPE* val,int* sort_order,idx_t nnz, int nmode)
{
	idx_t* tindex = (idx_t* ) malloc(nnz*nmode*sizeof(idx_t));
	TYPE* tval = (TYPE* ) malloc(nnz*sizeof(TYPE));
	for(int i=0; i< nnz; i++)
	{
		idx_t* nnz_ptr = tindex + i*nmode;
		int jump = (pindex[i] - index)/nmode;
		for(int j=0; j<nmode; j++)
		{
			int j_old = sort_order[j];
			nnz_ptr[j] = pindex[i][j_old];	 
		}
		tval[i] = val[jump];
	}

	for(int i=0; i< nnz*nmode; i++)
	{
		index[i] = tindex[i];
	}

	for(int i=0; i< nnz; i++)
	{
		val[i] = tval[i];
	}

	rem(tval);
	rem(tindex);
	return 0;
}

int read_tensor(const char* file, csf* res,  coo* debugt, int order_num)
{
	FILE *fp;
	//int *loc;
	#define  buf_size  (MAX_MODE+1)*30
	char buf[buf_size];
	//int nflag, sflag;
	// int dummy, pre_count=0, tmp_ne;
	int i,ii=-1;
	idx_t* index = NULL;
	idx_t** pindex = NULL;
	int sizestep = 1000000;
	int size = sizestep;
	int nmode = 0;
	int nnz = 0;
	int* mlen=NULL;
	idx_t * fiber_count;
	TYPE* vals=NULL;
	srand(time(NULL));
	//int compare_nnz(const void *a, const void *b);


	fp = fopen(file, "r");

	//for(int i = 0 ; i < 10 ; i++)
	while(fgets(buf, buf_size, fp))
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
		//printf("*%s* val is %lf mode len is %d\n",buf, val,mode_len);

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
	fiber_count = (idx_t*) malloc(nmode*sizeof(idx_t));

	if(order_num == -1)  // For SpTL1
		order_modes(mlen, nmode, sort_order);
	else if (order_num == -2) // for SpTL1so
	{
		order_modes(mlen, nmode, sort_order);
		int temp_mode_id = sort_order[nmode-1];
		sort_order[nmode-1] = sort_order[nmode-2];
		sort_order[nmode-2] = temp_mode_id;
		int temp_mode_len = mlen[nmode-1];
		mlen[nmode-1] = mlen[nmode-2];
		mlen[nmode-2] = temp_mode_len;
	}
	else if (order_num == -3) // for SpTL2
	{
		order_modes(mlen, nmode, sort_order);
		int temp_mode_id = sort_order[nmode-1];
		for (int i = nmode -2 ; i >= 0 ; i--)
		{
			sort_order[i+1] = sort_order[i];
		}
		sort_order[0] = temp_mode_id;
		int temp_mode_len = mlen[nmode-1];
		for (int i = nmode -2 ; i >= 0 ; i--)
		{
			mlen[i+1] = mlen[i];
		}
		mlen[0] = temp_mode_len;
	}
	else if (order_num == -4) // for SpTL2so
	{
		order_modes(mlen, nmode, sort_order);
		int temp_mode_id = sort_order[nmode-1];
		sort_order[nmode-1] = sort_order[nmode-2];
		sort_order[nmode-2] = temp_mode_id;
		int temp_mode_len = mlen[nmode-1];
		mlen[nmode-1] = mlen[nmode-2];
		mlen[nmode-2] = temp_mode_len;
		temp_mode_id = sort_order[nmode-1];
		for (int i = nmode -2 ; i >= 0 ; i--)
		{
			sort_order[i+1] = sort_order[i];
		}
		temp_mode_len = mlen[nmode-1];
		for (int i = nmode -2 ; i >= 0 ; i--)
		{
			mlen[i+1] = mlen[i];
		}
		mlen[0] = temp_mode_len;

	}
	else
	{
		idx_t* len = (idx_t*) malloc(sizeof(idx_t)*nmode);
		create_perm(order_num, sort_order, nmode);
		printf("Trying order ");
		for(int i = 0; i<nmode ; i++)
		{
			printf("-> %d ",sort_order[i] );
		}		
		printf("\n");

		for(int i = 0; i<nmode ; i++)
		{
			len[i] = mlen[sort_order[i]];
		}

		for(int i = 0; i<nmode ; i++)
		{
			mlen[i] = len[i];
		}
		rem(len);
	}

	printf("Tensor has dimensions ");
	for(i = 0; i<nmode-1 ; i++)
		printf("%dx",mlen[i]);
	printf("%d and %d nnz\n",mlen[i],nnz);
	/*
	sort_order[0] = 2;
	sort_order[1] = 0;
	sort_order[2] = 1;
	sort_order[3] = 3;
	*/
	// COO reading is finished
	{	
		auto start = std::chrono::high_resolution_clock::now();
		qsort(pindex,nnz,sizeof(idx_t*),compare_nnz);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> dif = end-start;
		
		printf("time for sorting COO is %lf \n",dif.count());
	}
	count_fiber(pindex,nnz,nmode,ii,fiber_count,sort_order);

	if(VERBOSE > VERBOSE_HIGH)
		print_stats(sort_order, fiber_count, file, nmode, nnz);


	//reorder_stat(nmode, nnz,  pindex, file );

	csf *t = res;


	coo2csf(pindex, index, vals, nnz,nmode,fiber_count,t,mlen,sort_order);
	if(VERBOSE == VERBOSE_DEBUG)
	{
		printt(pindex , nnz , nmode, vals);

		printf("nmode %d \n ",t->nmode);
	}


	if(debugt == NULL)
	{
		rem(index);
		
		rem(vals);
		
	}
	else
	{
		sort_coo(pindex,index,vals,sort_order,nnz,nmode);


		debugt -> ind = index;
		debugt -> val = vals;
		debugt -> nnz = nnz;
		debugt -> nmode = nmode;
		//debugt -> sort_order = sort_order;
	}
	rem(pindex);
	rem(sort_order);
	rem(mlen)
	rem(fiber_count);

	return 0;
}



#endif
