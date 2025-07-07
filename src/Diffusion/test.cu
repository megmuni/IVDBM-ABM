
#include "test.cuh"

__global__
void hello_kernel(char *a, int *b)
{
 a[threadIdx.x] += b[threadIdx.x];
}

void hello(char *ad, int *bd)
{
 dim3 dimBlock( blocksize, 1 );
 dim3 dimGrid( 1, 1 );
 hello_kernel<<<dimGrid, dimBlock>>>(ad, bd);
}

