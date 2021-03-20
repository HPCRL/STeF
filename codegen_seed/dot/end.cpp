// write to output matrix
TYPE* matval = mats[0]->val + ((mats[0]) -> dim2) * t->ind[0][i0]; 

#pragma omp simd
for(int y=0 ; y<r ; y++)
{
//	printf("0th level loop %lf\n",partial_results[y]);
matval[y] = partial_results[y];
partial_results[y] = 0;
}

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
rem(partial_results_all);	
return 0;

}