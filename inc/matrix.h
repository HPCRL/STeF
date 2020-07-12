#ifndef MATRIX_H
#define MATRIX_H
#include "../inc/util.h"
#define MAT(m,x,y) m->val[x*m->dim2 + y]

struct dense_matrix
{
	TYPE* val;
	idx_t dim1;
	idx_t dim2;
};

typedef struct dense_matrix matrix;

matrix* create_matrix(int dim1, int dim2, TYPE val = 0);

int print_matrix(matrix mat);
int free_matrix(matrix* mat);

#endif