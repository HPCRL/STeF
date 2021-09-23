#ifndef READER_H
#define READER_H
#include "../inc/util.h"
#include "../inc/tensor.h"


int read_tensor(const char* file,  csf* res = NULL, coo* debugt = NULL, int order_num = -1);
int coo2csf_tensor(const char* file, csf* res,  coo* debugt, int order_num);


#endif