#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h> 
#include <string.h>
#include <chrono>

//#define OMP

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


#ifdef OMP
#include <omp.h>
#endif

#define idx_t int
#define MAX_MODE 10
#define TYPE float



#define VERBOSE VERBOSE_LOW


enum verbosity
{
	VERBOSE_SILENT,
	VERBOSE_LOW,
	VERBOSE_HIGH,
	VERBOSE_DEBUG
};


#endif