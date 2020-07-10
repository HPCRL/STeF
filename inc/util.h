#ifndef UTIL_H
#define UTIL_H
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h> 
#include <string.h>
#include <chrono>

//#define OMP
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