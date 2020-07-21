#ifndef MTTKRP_HARDWIRED
#define MTTKRP_HARDWIRED
#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"

int mttkrp_hardwired_first_3(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_first_4(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_first_5(csf* t, int mode, int r, matrix** mats, int profile = -1);


#endif