#ifndef IVDBM_CHEMICAL_CHANNEL_VIEWS_H
#define IVDBM_CHEMICAL_CHANNEL_VIEWS_H

/**
 * @file chemical_channel_views.h
 * @brief Canonical names for per-species patch buffers (Phase III).
 *
 * During transition, secretion_delta and the diffusion export share the legacy
 * d* channel within one tick (diffuse then cells, then merge). Future work may
 * split into separate arrays.
 */

#include <memory>
#include <vector>

/** Owns one species' per-patch fields (Phase III+); optional shared storage. */
using ChemicalFieldBuffer = std::shared_ptr<std::vector<float>>;

/** Non-owning views into legacy chemAllocation rows for one species. */
struct SpeciesChannelViews
{
    float *concentration = nullptr;    /**< Legacy p* */
    float *secretion_delta = nullptr;  /**< Legacy d* (accumulator pre-merge) */
    float *diffused = nullptr;         /**< Same as secretion_delta until channels split */
};

#endif
