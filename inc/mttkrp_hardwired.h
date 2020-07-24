#ifndef MTTKRP_HARDWIRED
#define MTTKRP_HARDWIRED
#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"
#include "../inc/mutex.h"

int reduce(csf* t, int r, matrix* mat);
int mttkrp_hardwired_first_3(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_first_4(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_first_5(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_last_3(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile = -1);
int mttkrp_hardwired_last_4(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile = -1);
int mttkrp_hardwired_last_5(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile = -1);
int mttkrp_hardwired_last_vec_2(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile);
int mttkrp_hardwired_last_vec_3(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile);
int mttkrp_hardwired_last_vec_4(csf* t, int mode, int r, matrix** mats, mutex_array* mutex, int profile);


#endif