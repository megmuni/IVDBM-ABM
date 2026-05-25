#ifndef IVDBM_DIFFUSION_ALGORITHM_NAME_H
#define IVDBM_DIFFUSION_ALGORITHM_NAME_H

#include "../Diffusion3D/core/diffusion_algorithm.h"

/** Human-readable label for stdout / logs (matches Diffusion3D demo naming). */
const char *diffusion_algorithm_label(DiffusionAlgorithm algo);

/**
 * @brief Default patch-field solver: FFT when built with DIFFUSION3D_CUDA, else CPU explicit.
 */
DiffusionAlgorithm patch_field_default_diffusion_algorithm();

#endif
