// MemoryMappedAdapter: zero-copy float** <-> std::vector<double>
#ifndef DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H
#define DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H

/**
 * @file memory_mapped_adapter.h
 * @brief Zero-copy memory mapping utilities for Diffusion3D.
 *
 * Provides MemoryMappedAdapter for mapping between std::vector<double> and float* buffers
 * without copying, for efficient interop between legacy and modern code.
 */

#pragma once
#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

class MemoryMappedAdapter
{
public:
    // Map a std::vector<double> to a float** (no copy, reinterpret)
    static float *map_double_vector_to_float(std::vector<double> &data)
    {
        if (data.empty())
            throw std::runtime_error("Empty vector");
        return reinterpret_cast<float *>(data.data());
    }
    // Validate alignment for zero-copy
    static void validate_alignment(const std::vector<double> &data)
    {
        if (reinterpret_cast<uintptr_t>(data.data()) % alignof(float) != 0)
            throw std::runtime_error("Data not properly aligned for float mapping");
    }
};

#endif // DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H
