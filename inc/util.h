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
//#include <typeinfo>

//#define OMP



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

#define idx_t int
#define MAX_MODE 10
#define TYPE double
#define VERBOSE VERBOSE_HIGH
#define CORRECTNESS_THRESHOLD 1E-10
#define DOT_PARALLEL_DEPTH 2	
#define PAD 64

enum verbosity
{
	VERBOSE_SILENT,
	VERBOSE_LOW,
	VERBOSE_HIGH,
	VERBOSE_DEBUG
};



#endif
