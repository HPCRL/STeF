for(idx_t iYYY = t->ptr[XXX][iXXX] ; iYYY< t->ptr[XXX][iXXX+1]; iYYY++)
{
TYPE* pr = partial_results + XXX * r;
TYPE tval = t->val[iYYY];
TYPE* matval = (mats[YYY]->val) + ((mats[YYY]) -> dim2) * t->ind[YYY][iYYY];

#pragma omp simd
for(int y=0 ; y<r ; y++)
{
pr[y] += tval * matval[y];	// TTM step			
}
}	