#ifndef IVDBM_DECAY_CONSTANTS_H
#define IVDBM_DECAY_CONSTANTS_H

#include <cmath>

/**
 * @file decay_constants.h
 * @brief Clearance terms for one species over one tick. Defaults are inert.
 *
 */

struct DecayConstants {
    double k = 0.0;          /**< 1/min, = ln2 / half-life */
    double retain = 1.0;     /**< exp(-k*dt) applied to the stored field */
    double src_factor = 0.0; /**< (1-retain)/k in minutes; used once channels split */

    /** Build from a half-life in minutes. <= 0 yields inert constants */
    static DecayConstants from_half_life(double half_life_min, double dt_min) {
        DecayConstants d;
        if (half_life_min <= 0.0) {
            d.src_factor = dt_min;
            return d;
        }
        d.k = M_LN2 / half_life_min;
        d.retain = std::exp(-d.k * dt_min);
        d.src_factor = (1.0 - d.retain) / d.k;
        return d;
    }
};

#endif
