#include "../inc/reader.h"
#include "../inc/matrix.h"
#include "../inc/mttkrp.h"
#include "../inc/mttkrp_combined.h"
#include "../inc/mttkrp_hardwired.h"
#include <time.h>
#include <string>

int main(int argc, char** argv)
{
    mutex_array* mutex = mutex_alloc_custom(10000 , 16);

    double atomic_threshold = 0;
    int r = 32;
    if(argc > 1)
        r= atoi(argv[1]);
    atomic_threshold = atomic_thresh(r,mutex);
    printf("Atomic threshold is %lf\n",atomic_threshold);
    return 0;
}