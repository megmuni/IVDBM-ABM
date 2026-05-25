#include "diffusion_algorithm_name.h"

const char *diffusion_algorithm_label(DiffusionAlgorithm algo)
{
    switch (algo)
    {
    case DiffusionAlgorithm::ExplicitHeatEquation:
        return "ExplicitHeatEquation (CPU stencil)";
    case DiffusionAlgorithm::GpuStencil:
        return "GpuStencil (GPU stencil)";
    case DiffusionAlgorithm::GpuFftPrecomputed:
        return "GpuFftPrecomputed (FFT)";
    }
    return "unknown diffusion algorithm";
}

DiffusionAlgorithm patch_field_default_diffusion_algorithm()
{
#ifdef DIFFUSION3D_CUDA
    return DiffusionAlgorithm::GpuFftPrecomputed;
#else
    return DiffusionAlgorithm::ExplicitHeatEquation;
#endif
}
