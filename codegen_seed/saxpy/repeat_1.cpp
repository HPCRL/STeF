for(idx_t iYYY = t->ptr[XXX][iXXX]; iYYY < t->ptr[XXX][iXXX+1] ; iYYY++)	
{
TYPE* matvalYYY = mats[YYY]->val + ((mats[YYY]->dim2) * t->ind[YYY][iYYY]);
//printf(" middle index is %d\n",i1);

#pragma omp simd
for(int y=0; y<r ; y++)
{
partial_products[y+YYY*r] = partial_products[y+XXX*r] * matvalYYY[y];	
}