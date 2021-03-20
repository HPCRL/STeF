int mttkrp_hardwired_first_not_fused_NNN(csf* t, int mode, int r, matrix** mats, int profile )
{
int nmode = t->nmode;
int num_th = 1;
int partial_results_size = nmode*r+PAD;
#ifdef OMP
num_th = omp_get_max_threads();

#endif

printf("num ths %d\n", num_th);

TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

for(int i = 0 ; i < num_th*partial_results_size ; i++)
partial_results_all[i] = 0;

set_matrix(*mats[0],0);

#ifdef OMP
#pragma omp parallel 
#endif
{
int th = 0;
#ifdef OMP
th = omp_get_thread_num();
#endif
auto time_start = std::chrono::high_resolution_clock::now();
TYPE* partial_results = partial_results_all + th*partial_results_size;
#ifdef OMP
#pragma omp for schedule(dynamic,1)
#endif
for(idx_t i0 = 0 ; i0< t->fiber_count[0]; i0++)
{