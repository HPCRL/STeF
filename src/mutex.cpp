#ifndef MUTEX_CPP
#define MUTEX_CPP
#include "../inc/util.h"
#ifdef OMP
#include "../inc/mutex.h"

mutex_array * mutex_alloc_custom(int const num_locks,int const pad_size)
{
	mutex_array * array = (mutex_array *) malloc(sizeof(mutex_array));

	array->num_locks = num_locks;
	array->pad_size = pad_size;


	array->locks = (omp_lock_t *) malloc(num_locks * pad_size * sizeof(*array->locks));
	for(int l=0; l < num_locks; ++l) 
	{
		int const lock = mutex_translate_id(l, num_locks, pad_size);
		omp_init_lock(array->locks + lock);
	}

	return array;
}
#endif
#endif
