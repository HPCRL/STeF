#ifndef MTTKRP_H
#define MTTKRP_H
#include "../inc/matrix.h"
#include "../inc/tensor.h"


int mttkrp_atomic3(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic(csf* t,int mode, int r, matrix** mats);


#endif