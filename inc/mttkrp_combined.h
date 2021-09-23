#ifndef MTTKRP_COMBINED
#define MTTKRP_COMBINED

#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"
#include "../inc/mutex.h"
#include "../inc/mttkrp_hardwired.h"


template <int mode, bool intv1>
int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile );

int mttkrp_combined_3(csf* t, int r, matrix** mats, int profile , int mode, bool intv1 );

template <int mode, bool intv1, bool intv2>
int mttkrp_combined_4(csf* t, int r, matrix** mats, int profile );

template <int mode, bool intv1, bool intv2, bool intv3>
int mttkrp_combined_5(csf* t, int r, matrix** mats, int profile );

#endif