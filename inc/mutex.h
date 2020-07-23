#ifndef MUTEX_H
#define MUTEX_H
#include "../inc/util.h"

#ifdef OMP
#include <omp.h>
#endif
typedef struct
{
  bool initialized;
  int num_locks;
  int pad_size;
  #ifdef OMP
  omp_lock_t * locks;
  #else
  int* locks;
  #endif
} mutex_array;


mutex_array * mutex_alloc_custom(int const num_locks,int const pad_size);

void mutex_free(mutex_array * array);



static inline int mutex_translate_id(int const id, int const num_locks,int const pad_size)
{
  return (id & 0x3FF)>>4; //(id % num_locks) * pad_size;
}


static inline void mutex_set_lock(mutex_array * const array,int const id)
{
  int const lock_id = mutex_translate_id(id, array->num_locks, array->pad_size);
  #ifdef OMP
  omp_set_lock(array->locks + lock_id);
  #endif
}

static inline void mutex_unset_lock(mutex_array * const array,int const id)
{
  int const lock_id = mutex_translate_id(id, array->num_locks, array->pad_size);
  #ifdef OMP
  omp_unset_lock(array->locks + lock_id);
  #endif
}






#endif