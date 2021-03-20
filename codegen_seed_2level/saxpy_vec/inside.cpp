for(idx_t iYYY = t->ptr[XXX][iXXX]; iYYY < t->ptr[XXX][iXXX+1]  ; iYYY++)
{
const idx_t row_id = t->ind[YYY][iYYY];
TYPE* xx  = vals + t->ind[YYY][iYYY]*(mats[mode]->dim2);
TYPE* yy = partial_products + XXX*r ;
TYPE* intval = t->intval[YYY] + iYYY*r ;


#pragma omp simd
for(int i=0 ; i<r ; i++)
{
TYPE increment = yy [i] * intval[i];
xx [i]	+= increment;
}

}