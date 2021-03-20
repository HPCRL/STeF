}
auto time_end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> time_diff = time_end-time_start;
printf("Hardwired time for mode %d thread %d %lf \n",t->modeid[mode],th,time_diff.count());		

if(profile == mode)
{
LIKWID_MARKER_STOP("Compute");
}
}

LIKWID_MARKER_CLOSE;

t->num_th = num_th;
if(num_th > 1)
{
auto time_start = std::chrono::high_resolution_clock::now();
reduce(t,r,mats[mode]);
auto time_end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> time_diff = time_end-time_start;
printf("Hardwired time for reducing mode %d is %lf \n",t->modeid[mode],time_diff.count());		
}

rem(partial_products_all);
return 0;
}