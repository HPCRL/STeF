#ifndef READER_H
#define READER_H
#include "../inc/util.h"
#include "../inc/tensor.h"


int read_tensor(const char* file,  csf* res = NULL, coo* debugt = NULL, int order_num = -1);

#endif