#include <cassert>

#include "diffusion3d_gpu.h"

void diffusion3d_step_euler_gpu(int, int, int, double, double, double, double *&, double *&,
                                int, int, int)
{
    assert(false && "diffusion3d_step_euler_gpu requires DIFFUSION3D_CUDA build");
}
