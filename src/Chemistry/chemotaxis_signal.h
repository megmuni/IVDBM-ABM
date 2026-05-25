#ifndef IVDBM_CHEMOTAXIS_SIGNAL_H
#define IVDBM_CHEMOTAXIS_SIGNAL_H

/**
 * @file chemotaxis_signal.h
 * @brief Per-patch signal used for cell movement toward higher cytokine levels.
 */

/** Non-owning view of chemotaxis strength on each patch (length nx×ny×nz). */
struct ChemotaxisSignal
{
    float *data = nullptr;
};

#endif
