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
		//random_matrix(*mats[i],i);
		set_matrix(*mats[i],1);
		if(VERBOSE  == VERBOSE_DEBUG)
			print_matrix(*mats[i]);
	}

	mttkrp_fused_init(t,r,true);
	b_thread_start(t);

	double total=0;

	int num_th = 1;

	#ifdef OMP
	num_th = omp_get_max_threads();
	#endif

	int run_saved = 0;
	if (argc > 5)
		run_saved = atoi(argv[5]);
		
	if (nmode == 3)
	{
		if(run_saved != 2)
		for(int mode = 0 ; mode<nmode ; mode++)
		{
			const bool intv = false;
			clock_t cstart, cend;
			bool is_atomic = mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) );
			auto start = std::chrono::high_resolution_clock::now();
			cstart = clock();
			/*
			if(is_atomic)	
			{
				if (mode == 0)
					mttkrp_combined_3<0,false,false>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_3<1,false,false>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_3<2,false,false>(t,r,mats,profile);
			}
			else
			{
				if (mode == 0)
					mttkrp_combined_3<0,false,true>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_3<1,false,true>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_3<2,false,true>(t,r,mats,profile);
			}
			*/	
			const bool saved = true;

			if(is_atomic)	
			{
				if (mode == 0)
					mttkrp_combined_lb_3<0,saved,false>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_lb_3<1,saved,false>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_lb_3<2,saved,false>(t,r,mats,profile);
			}
			else
			{
				if (mode == 0)
					mttkrp_combined_lb_3<0,saved,true>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_lb_3<1,saved,true>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_lb_3<2,saved,true>(t,r,mats,profile);
			}
			cend = clock();
			//printf("here\n");
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			total += diff.count();
			printf("Combined %s not saved time for mode %d %lf \n",(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
			double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
			// printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
			if(debug)
			{
				/*
				set_matrix(*mats[0],1);
				set_matrix(*mats[1],1);
				set_matrix(*mats[2],1);
				*/
				auto start2 = std::chrono::high_resolution_clock::now();
				mttkrp_test(dt,mode,r,mats);
				auto end2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end2-start2;
				//total += diff.count();
				// printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}

			for(i=0 ; i<nmode ; i++)
			{
				
				if(VERBOSE  == VERBOSE_DEBUG)
				{
					print_matrix(*mats[i]);
					
				}
				//random_matrix(*mats[i],i+1);
				//set_matrix(*mats[i],1);
			}
			//random_matrix(*mats[mode],mode);
			set_matrix(*mats[mode],1);
			set_matrix(*mats[0],1);
			set_matrix(*mats[1],1);
			set_matrix(*mats[2],1);

			
		}
		printf("Total Intermediate %s template MTTKRP time %lf\n",("not saved"),total);
		memset(t->intval[1],0,sizeof(TYPE)*t->fiber_count[1]);
		total = 0;
		if(run_saved != 1)
		for(int mode = 0 ; mode<nmode ; mode++)
		{
			const bool intv = true;
			bool is_atomic = mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) );
			clock_t cstart, cend;
			auto start = std::chrono::high_resolution_clock::now();
			cstart = clock();
			if(is_atomic)	
			{
				
				if (mode == 0)
					mttkrp_combined_3<0,true,false>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_3<1,true,false>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_3<2,true,false>(t,r,mats,profile);
			}
			else
			{
				if (mode == 0)
					mttkrp_combined_3<0,true,true>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_3<1,true,true>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_3<2,true,true>(t,r,mats,profile);
			}	
			cend = clock();
			//printf("here\n");
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			total += diff.count();
			printf("Combined %s saved time for mode %d %lf \n",(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
			double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
			// printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
			if(debug)
			{
				auto start2 = std::chrono::high_resolution_clock::now();
				mttkrp_test(dt,mode,r,mats);
				auto end2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end2-start2;
				//total += diff.count();
				// printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}

			for(i=0 ; i<nmode ; i++)
			{
				
				if(VERBOSE  == VERBOSE_DEBUG)
				{
					print_matrix(*mats[i]);
					
				}
				//random_matrix(*mats[i],i+1);
				//set_matrix(*mats[i],1);
			}
			random_matrix(*mats[mode],mode);
			
		}
		printf("Total Intermediate %s template MTTKRP time %lf\n",( "saved" ),total);
	}
	else if (nmode == 4)
	{
		int num_cases = 1;
		for(int i = 0 ; i< 4-2; i++)
			num_cases *= 2;
		total = 0;
		for(int mode = 0 ; mode<nmode ; mode++)
		{
			const bool intv1 = false;
			const bool intv2 = true;
			const int mode_c = mode;
			bool is_atomic = mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) );
			clock_t cstart, cend;
			auto start = std::chrono::high_resolution_clock::now();
			cstart = clock();
			if(is_atomic)	
			{	// Atomic Update
				if (mode == 0)
					mttkrp_combined_4<0,intv1,intv2,false>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_4<1,intv1,intv2,false>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_4<2,intv1,intv2,false>(t,r,mats,profile);
				else if (mode == 3)
					mttkrp_combined_4<3,intv1,intv2,false>(t,r,mats,profile);
			}
			else
			{ 	// Privatization
				if (mode == 0)
					mttkrp_combined_4<0,intv1,intv2,true>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_4<1,intv1,intv2,true>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_4<2,intv1,intv2,true>(t,r,mats,profile);
				else if (mode == 3)
					mttkrp_combined_4<3,intv1,intv2,true>(t,r,mats,profile);
			}	
			cend = clock();
			//printf("here\n");
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			total += diff.count();
			printf("Combined %s not saved time for mode %d %lf \n",(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
			double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
			//printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
			if(debug)
			{
				auto start2 = std::chrono::high_resolution_clock::now();
				mttkrp_test(dt,mode,r,mats);
				auto end2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end2-start2;
				//total += diff.count();
				// printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}

			for(i=0 ; i<nmode ; i++)
			{
				
				if(VERBOSE  == VERBOSE_DEBUG)
				{
					print_matrix(*mats[i]);
					
				}
				//random_matrix(*mats[i],i+1);
				//set_matrix(*mats[i],1);
			}
			random_matrix(*mats[mode],mode);

			
		}

		printf("Total Intermediate %s template MTTKRP time %lf\n",("not saved"),total);
	}
	else if (nmode == 5)
	{
		int num_cases = 1;
		for(int i = 0 ; i< 5-2; i++)
			num_cases *= 2;
		total = 0;
		for(int mode = 0 ; mode<nmode ; mode++)
		{
			const bool intv1 = false;
			const bool intv2 = false;
			const bool intv3 = false;
			const int mode_c = mode;
			bool is_atomic = mode > 0 && ((t->fiber_count[mode] / t->mlen[mode] < num_th * ATOMIC_THRESH) || (t->mlen[mode] * r * num_th >= PRIVATIZED_THRESH) );
			clock_t cstart, cend;
			auto start = std::chrono::high_resolution_clock::now();
			cstart = clock();
			if(is_atomic)	
			{	// Atomic Update
				const bool privatized = false;
				if (mode == 0)
					mttkrp_combined_5<0,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_5<1,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_5<2,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 3)
					mttkrp_combined_5<3,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 4)
					mttkrp_combined_5<4,intv1,intv2,intv3,privatized>(t,r,mats,profile);
			}
			else
			{ 	// Privatization
				const bool privatized = true;
				if (mode == 0)
					mttkrp_combined_5<0,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 1)
					mttkrp_combined_5<1,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 2)
					mttkrp_combined_5<2,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 3)
					mttkrp_combined_5<3,intv1,intv2,intv3,privatized>(t,r,mats,profile);
				else if (mode == 4)
					mttkrp_combined_5<4,intv1,intv2,intv3,privatized>(t,r,mats,profile);
			}	
			cend = clock();
			//printf("here\n");
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			total += diff.count();
			printf("Combined %s not saved time for mode %d %lf \n",(is_atomic ? "atomic    " : "privatized"),t->modeid[mode],diff.count());
			double cdiff = ((double) (cend - cstart)) / CLOCKS_PER_SEC;
			//printf("Clock time for mode %d is %lf \n",t->modeid[mode],cdiff);
			if(debug)
			{
				auto start2 = std::chrono::high_resolution_clock::now();
				mttkrp_test(dt,mode,r,mats);
				auto end2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end2-start2;
				//total += diff.count();
				// printf("COO sequential time for mode %d %lf \n",t->modeid[mode],diff.count());
			}

			for(i=0 ; i<nmode ; i++)
			{
				
				if(VERBOSE  == VERBOSE_DEBUG)
				{
					print_matrix(*mats[i]);
					
				}
				//random_matrix(*mats[i],i+1);
				//set_matrix(*mats[i],1);
			}
			random_matrix(*mats[mode],mode);

			
		}

		printf("Total Intermediate %s template MTTKRP time %lf\n",("not saved"),total);
	}
	

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
