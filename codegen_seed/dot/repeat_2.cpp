// write to intval
TYPE* matval = (mats[YYY]->val) + ((mats[YYY]) -> dim2) * t->ind[YYY][iYYY];
TYPE* intval = t->intval[YYY] + iYYY*r;

#pragma omp simd
for(int y=0 ; y<r ; y++)
{
partial_results[y + XXX*r] += partial_results[YYY*r+y] * matval[y]; // TTV
intval[y] = partial_results[YYY*r + y];
}

#pragma omp simd
for(int y=0; y<r; y++)
{
partial_results[YYY*r+y] = 0;
}
}