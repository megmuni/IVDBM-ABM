#ifndef DIFFUSION3D_DIFFUSION_ALGORITHM_H
#define DIFFUSION3D_DIFFUSION_ALGORITHM_H

/**
 * @file diffusion_algorithm.h
 * @brief PDE stepping algorithm selection for MultiSpeciesDiffusionEngine.
 */

enum class DiffusionAlgorithm
{
    ExplicitHeatEquation = 0,
    GpuStencil = 1,
    GpuFftPrecomputed = 2,
};

#endif // DIFFUSION3D_DIFFUSION_ALGORITHM_H
