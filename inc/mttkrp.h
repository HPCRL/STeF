#ifndef MTTKRP_H
#define MTTKRP_H
#include "../inc/matrix.h"
#include "../inc/tensor.h"
#include "../inc/util.h"
#include "../inc/mttkrp_hardwired.h"
#include "../inc/mutex.h"



int mttkrp_atomic3(csf* t, int mode, int r, matrix** mats);
int mttkrp_atomic_last(csf* t, int mode, int r, matrix** mats, int vec , int profile = -1);
int mttkrp_private_last(csf* t, int mode, int r, matrix** mats, int vec, int profile = -1);
int mttkrp_atomic_last_vec(csf* t, int mode, int r, matrix** mats, int vec , int profile = -1);
int mttkrp_private_last_vec(csf* t, int mode, int r, matrix** mats, int vec, int profile = -1);
int mttkrp_atomic_first(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_atomic_middle(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_atomic(csf* t,int mode, int r, matrix** mats, int profile = -1);
int mttkrp_private(csf* t,int mode, int r, matrix** mats, int profile = -1);
int mttkrp_fused_init(csf* t,int r, bool cap=false);
int mttkrp_fused_free(csf* t,int r);
int mttkrp_test(coo* dt, int mode, int r, matrix** mats);
int find_inds(idx_t* inds ,csf* t,idx_t it);
int dist_dot_work(idx_t* inds ,csf* t,int p,idx_t* count,int th,int depth=DOT_PARALLEL_DEPTH);
int mttkrp_hardwired_first_not_fused(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_first(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired_last(csf* t, int mode, int r, matrix** mats, int profile = -1);
int mttkrp_hardwired(csf* t, int mode, int r, matrix** mats, int profile);
int mttkrp_hardwired_middle(csf* t, int mode, int r, matrix** mats, int profile);
#endif