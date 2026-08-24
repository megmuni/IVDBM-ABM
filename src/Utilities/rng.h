#ifndef ABM_UTILITIES_RNG_H
#define ABM_UTILITIES_RNG_H

/**
 * rng.h
 *
 * File contents: the only source of randomness in the model. Every random
 * decision draws from a Philox stream keyed by the world seed and a stable
 * per-entity id, so a run reproduces across machines and thread counts.
 *
 * Philox is a counter-based generator: a draw is a pure function of
 * (seed, entity id, draw index). Original source:
 * https://doi.org/10.1145/2063384.2063405
 */

#include <cstdint>
#include <iterator>

#include <Random123/philox.h>

namespace abm {
namespace rng {

/** The world seed. Set once at startup, before any Agent or ECM is built. */
inline uint32_t &global_seed() {
  static uint32_t seed = 0;
  return seed;
}

inline void set_global_seed(uint32_t seed) { global_seed() = seed; }

/**
 * Stream: one entity's random sequence. Stored by value; holds no pointers and
 * shares no state, so agents on different threads never interact through it.
 *
 * The two ids identify the entity. They must be derived from deterministic
 * state (a patch index, a birth tick), never from a counter handed out in
 * thread arrival order.
 */
class Stream {
public:
  typedef uint32_t result_type;

  Stream() : substate_(4) {
    counter_[0] = 0;
    counter_[1] = 0;
    counter_[2] = 0;
    counter_[3] = 0;
    key_[0] = 0;
    key_[1] = 0;
  }

  Stream(uint32_t id0, uint32_t id1) : substate_(4) {
    counter_[0] = 0;
    counter_[1] = id0;
    counter_[2] = id1;
    counter_[3] = 0;
    key_[0] = global_seed();
    key_[1] = 0;
  }

  static result_type min() { return 0; }
  static result_type max() { return 0xFFFFFFFFu; }

  result_type operator()() {
    if (substate_ >= 4) {
      result_ = r123::Philox4x32()(counter_, key_);
      ++counter_[0];
      substate_ = 0;
    }
    return result_[substate_++];
  }

private:
  r123::Philox4x32::ctr_type counter_;
  r123::Philox4x32::ukey_type key_;
  int substate_;
  r123::Philox4x32::ctr_type result_;
};

/** A uniform draw in [0, 1). */
inline double uniform(Stream &stream) {
  return stream() * (1.0 / 4294967296.0);
}

/** A uniform integer in [0, n). */
inline int uniform_int(Stream &stream, int n) {
  return static_cast<int>(stream() % static_cast<uint32_t>(n));
}

/** True with the given probability, expressed as a percentage in [0, 100]. */
inline bool roll_percent(Stream &stream, double percent) {
  return uniform(stream) * 100.0 < percent;
}

/**
 * Fisher-Yates over [first, last). Written out rather than calling
 * std::shuffle, whose mapping from generator output to permutation differs
 * between standard library implementations and would break reproducibility
 * across compilers.
 */
template <typename RandomIt>
void shuffle(Stream &stream, RandomIt first, RandomIt last) {
  typename std::iterator_traits<RandomIt>::difference_type n = last - first;
  for (typename std::iterator_traits<RandomIt>::difference_type i = n - 1; i > 0;
       --i) {
    typename std::iterator_traits<RandomIt>::difference_type j =
        static_cast<typename std::iterator_traits<RandomIt>::difference_type>(
            stream() % static_cast<uint32_t>(i + 1));
    typename std::iterator_traits<RandomIt>::value_type tmp = first[i];
    first[i] = first[j];
    first[j] = tmp;
  }
}

} // namespace rng
} // namespace abm

#endif /* ABM_UTILITIES_RNG_H */
