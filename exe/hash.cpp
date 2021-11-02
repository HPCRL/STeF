#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"

inline idx_t  hash_index(idx_t* nnz, idx_t* pad, int* modes, int nmodes, idx_t prime)
{
    idx_t hash_val = 0;
    for (int i = 0 ; i<nmodes ; i++)
    {
        hash_val += nnz[modes[i]]*pad[i];
        //printf("%lld-%lld " ,nnz[modes[i]],pad[i]);
    }
    //printf("%lld\n",hash_val);
    return hash_val % prime;
}

idx_t hash_dim(coo* dt,int* modes, int nmodes)
{
    idx_t* pad = new idx_t[nmodes];
    pad[0] = 1;
    for (int i = 0 ; i<nmodes-1 ; i++)
    {
        pad[i+1] = pad[i]*dt->mlen[modes[i]];
    }

    idx_t prime = dt->nnz+1;

    idx_t* hash_table = new idx_t[dt->nnz];
    memset(hash_table,0,sizeof(idx_t)*dt->nnz);

    idx_t* nnz_ptr = dt->ind;
    idx_t pro_cnt = 0;
    #pragma omp parallel for reduction(+:pro_cnt)
    for (int i=0; i<dt->nnz; i++)
    {
        idx_t hv = hash_index(dt->ind + (i * dt->nmode),pad,modes,nmodes,prime);        
        if(hash_table[hv] == 0)
        {
            hash_table[hv] = 1;
            pro_cnt += 1;
        }            
    }

    return pro_cnt;
}

int main(int argc, char** argv)
{


	if(argc == 1)
	{
		printf("Usage is %s <matrix name> <number of ranks (optional)> \
<kernel (-1: SpTL , -2: SpTLso , 0-d!: permutation of modes) (optional)> \
<mode to profile (optional) >  \n", argv[0] );
		exit(1);
	}

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

   	print_csf(t,argv[1]);
    
    
    for(int nmodes = dt->nmode-1;nmodes > 0 ; nmodes --)
    {
        int* modes = new int[nmodes];
        int* pad = new int[nmodes];
        int num_configs = 1;
        for(int i=nmodes - 1 ; i>=0 ; i--)
        {   
            pad[i] = num_configs;         
            num_configs *= dt->nmode - i;            
        }

        for(int i = 0 ; i< num_configs ; i++)
        {                  
            int rem = i;          
            printf("projection order is " );
            for (int m = 0 ; m < nmodes ; m++)
            {
                modes[m] = rem / pad[m];
                //printf("\n??? is %d\n",modes[m]);
                if (m > 0)
                {
                    modes[m] = (modes[m] + modes[m-1]+1) % dt->nmode;
                    //printf("\n??? is %d\n",modes[m]);
                    bool cont = true;
                    while(cont)
                    {
                        cont = false;
                        for(int j = 0 ; j < m ; j++)
                        {   
                            if (modes[m] == modes[j])
                                cont = true;
                        }
                        if (cont)
                        {
                            modes[m] = (modes[m] + 1) % dt->nmode;
                        }
                    }
                }
                //printf("\nrem is %d\n",rem);
                printf("%d ",modes[m]);
                rem = rem % pad[m];
            }
            //printf("\n");
            auto start = std::chrono::high_resolution_clock::now();
            idx_t res = hash_dim(dt,modes,nmodes);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end-start;
            double time = diff.count();

            printf("Hash estimation took %lf seconds and result is %lld\n",time,res);
        }
        delete [] modes;
    }
	
	return 0;
}
