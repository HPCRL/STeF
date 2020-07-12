#ifndef MTTKRP_H
#define MTTKRP_H
#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"


int mttkrp_atomic3(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats, int vec = 0);
int mttkrp_atomic_first(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic_middle(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic(csf* t,int mode, int r, matrix** mats);
int mttkrp_fused_init(csf* t,int r);


#endif