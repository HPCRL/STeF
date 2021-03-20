int mttkrp_hardwired_last_vec_NNN(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile)
{
TYPE* partial_products_all;
int nmode;
int num_th;
nmode = t->nmode;

#ifdef OMP
num_th = omp_get_max_threads();
#else
num_th = 1;
int th = 0;
#endif

int partial_products_size = nmode*r + PAD;

printf("num ths %d\n", num_th);
partial_products_all = (TYPE* ) malloc(num_th*(partial_products_size)*sizeof(TYPE));

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
#ifdef OMP
int th = omp_get_thread_num();
if(VERBOSE == VERBOSE_HIGH)
printf("th id is %d\n",th);
#endif

TYPE* partial_products;	
partial_products = partial_products_all + th*partial_products_size;

TYPE* vals;
if(t->num_th > 1)
{
vals = t->private_mats[th]->val;
}
else
{
vals = mats[mode]->val;
}

memset(vals, 0 , mats[mode]->dim1*mats[mode]->dim2*sizeof(TYPE));

auto time_start = std::chrono::high_resolution_clock::now();
#ifdef OMP
#pragma omp for schedule(dynamic,1)
#endif
for(idx_t i0 = 0; i0 < (t->fiber_count[0]) ; i0++)
{

TYPE* matval0 = mats[0]->val + ((mats[0]->dim2) * (t->ind)[0][i0]);

#pragma omp simd
for(int y=0; y<r ; y++)
{
partial_products[y] = matval0[y];	
}