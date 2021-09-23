#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/mttkrp_combined.h"
#include "../inc/mttkrp_hardwired.h"
#include <time.h>

int main(int argc, char** argv)
{
	int nmode,i,r;
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

	t->intval = NULL;


	print_csf(t,argv[1]);
	
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


    //mttkrp_combined_3<0,true>(t,r,mats,profile);
    //mttkrp_combined_3(t,r,mats,profile,0,true);
    /*
	for(mode = 0 ; mode<nmode ; mode++)
	{
		clock_t cstart, cend;
		auto start = std::chrono::high_resolution_clock::now();
		cstart = clock();
        if (mode == 0)
		    mttkrp_combined_3<0,true>(t,r,mats,profile);
        else if (mode == 1)
            mttkrp_combined_3<1,true>(t,r,mats,profile);
        else if (mode == 2)
            mttkrp_combined_3<2,true>(t,r,mats,profile);
		cend = clock();
		//printf("here\n");
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		total += diff.count();
		printf("Hardwired time for mode %d %lf \n",t->modeid[mode],diff.count());
		double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
		printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
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


	}*/
    for (int repeat = 0 ; repeat < 1 ; repeat ++)
    {
        bool intv = false;
        if ( repeat)
            intv = true;
        for(int mode = 0 ; mode<nmode ; mode++)
        {
            clock_t cstart, cend;
            auto start = std::chrono::high_resolution_clock::now();
            cstart = clock();
            if (mode == 0)
                mttkrp_combined_3(t,r,mats,profile,0,intv);
            else if (mode == 1)
                mttkrp_combined_3(t,r,mats,profile,1,intv);
            else if (mode == 2)
                mttkrp_combined_3(t,r,mats,profile,2,intv);
            cend = clock();
            //printf("here\n");
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end-start;
            total += diff.count();
            printf("Intermediate save is %s, combined time for mode %d %lf \n",(intv ? "on" : "off"),t->modeid[mode],diff.count());
            double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
            printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
            if(debug)
            {
                //  for(i=0 ; i<nmode ; i++)  if (i != mode)               random_matrix(*mats[i],i);
                auto start2 = std::chrono::high_resolution_clock::now();
                mttkrp_test(dt,mode,r,mats);
                auto end2 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> diff = end2-start2;
                //total += diff.count();
                printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
            }
            random_matrix(*mats[mode],mode);

            for(i=0 ; i<nmode ; i++)
            {
                //random_matrix(*mats[i],i);
                
                if(VERBOSE  == VERBOSE_DEBUG)
                {
                    print_matrix(*mats[i]);
                    
                }
            }
        }
    }
	printf("Total Hardwired MTTKRP time %lf\n",total);
	

	
	

	//free(t->intval[1]);
	free_csf(t);
	free_coo(dt);
	//rem(t);
	//rem(dt);
	for(i=0 ; i<nmode ; i++)
	{
		free_matrix(mats[i]);
	}
	rem(mats);
	return 0;
}
