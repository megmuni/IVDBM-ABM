#ifndef IVDBM_CHEMICAL_CHANNEL_VIEWS_H
#define IVDBM_CHEMICAL_CHANNEL_VIEWS_H

/**
 * @file chemical_channel_views.h
 * @brief Named views into per-species patch chemistry buffers.
 *
 * Within one tick, the secretion_delta and diffused rows hold cell secretion and 
 * diffusion output respectively, before merge adds it to the concentration row.
 */

#include <memory>
#include <vector>

/** Optional owned buffer for a single species grid (future shared storage). */
using ChemicalFieldBuffer = std::shared_ptr<std::vector<float>>;

/** Non-owning pointers into ChemicalEnvironment channel rows for one species. */
struct SpeciesChannelViews
{
    float *concentration = nullptr;   /**< Stored level (p* channel). */
    float *secretion_delta = nullptr; /**< Cell secretion this tick (s* channel). */
    float *diffused = nullptr;        /**< Diffusion delta this tick (d* channel). */
};

#endif
