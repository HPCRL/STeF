int mttkrp_hardwired_first_NNN(csf* t, int mode, int r, matrix** mats, int profile )
{
int nmode = t->nmode;
int num_th = 1;
int partial_results_size = nmode*r+PAD;
#ifdef OMP
num_th = omp_get_max_threads();

#endif

//mutex_array* mutex = mutex_alloc_custom((t->mlen)[0] , 16);

printf("num ths %d\n", num_th);

TYPE* partial_results_all = (TYPE*) malloc(num_th*partial_results_size*sizeof(TYPE));

for(int i = 0 ; i < num_th*partial_results_size ; i++)
partial_results_all[i] = 0;

set_matrix(*mats[0],0);

if(profile == mode)
{
printf("profiling mode %d  == %d\n",profile, t->modeid[profile] );
LIKWID_MARKER_INIT;
}

#ifdef OMP
#pragma omp parallel
#endif
{
if (profile == mode)
{
LIKWID_MARKER_THREADINIT;	
}
}

#ifdef OMP
#pragma omp parallel
#endif
{
if(profile == mode)
{
LIKWID_MARKER_START("Compute");
}
}

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
for(idx_t i1 = 0 ; i1< t->fiber_count[1]; i1++)
{