#ifndef IVDBM_CHEMOTAXIS_SIGNAL_H
#define IVDBM_CHEMOTAXIS_SIGNAL_H

/**
 * @file chemotaxis_signal.h
 * @brief Agent chemotaxis field on the patch grid (legacy pcellgrad).
 */

/** Non-owning view of per-patch chemotaxis strength (size nx * ny * nz). */
struct ChemotaxisSignal
{
    float *data = nullptr;
};

#endif
