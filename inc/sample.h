#ifndef SAMPLE_H
#define SAMPLE_H 
#include "../inc/util.h"
#include "../inc/tensor.h"
#include <stdlib.h>

int estimate_fiber(coo* dt, int* order, int* fiber_count, double sample_rate= 0.001);

#endif