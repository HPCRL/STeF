#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h> 
#include <string.h>
#include <chrono>
#include <thread>         
#include <mutex>  
#include <cmath>

//#include <typeinfo>

#define OMP



#define rem(X) if(X != NULL){free(X); X = NULL;}

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

//#define OMP
#ifdef OMP
#include <omp.h>
#endif

#define idx_t long long int
#define MAX_MODE 10
#define TYPE double
#define VERBOSE VERBOSE_LOW
#define CORRECTNESS_THRESHOLD 1E-10
#define DOT_PARALLEL_DEPTH 1	
#define PAD 0
#define ATOMIC_THRESH 5
#define PRIVATIZED_THRESH 80000000000 // 80 GB
//#define PRIVATIZED_THRESH 4000000000 // 4GB

enum verbosity
{
	VERBOSE_SILENT,
	VERBOSE_LOW,
	VERBOSE_HIGH,
	VERBOSE_DEBUG
};



#endif
