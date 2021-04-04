#include "../inc/reader.h"
#include "../inc/tensor.h"

int main(int argc, char** argv)
{
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

	auto start = std::chrono::high_resolution_clock::now();
	int* sort_order = new int[t->nmode];
	for(int i=0; i< t->nmode ; i++)
		sort_order[i] = i;
	count_fiber(dt,sort_order,dt->nmode);

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end-start;
	
	printf("Hash based fiber counting took %lf seconds \n",diff.count());	


	return 0;
}