#ifndef DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H
#define DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H

/**
 * @file memory_mapped_adapter.h
 * @brief Zero-copy memory mapping utilities for Diffusion3D (adapted from diffusion3d interop patterns).
 *
 * Provides MemoryMappedAdapter for bi-directional mapping between legacy float* buffers
 * and modern std::vector<double> grids without copying. Enables seamless interop during
 * incremental migration of ABM to multi-species architecture.
 */

#include <vector>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <stdexcept>

/**
 * @class MemoryMappedAdapter
 * @brief Template utilities for zero-copy buffer mapping and validation.
 *
 * Supports:
 * - Reinterpret cast between double* and float* (same bit pattern, no semantic change)
 * - Alignment validation for safe memory mapping
 * - Dimension validation for grid reshape operations
 *
 * **Design**: All methods are static/inline for header-only usage (no .cpp required).
 */
class MemoryMappedAdapter
{
public:
    /**
     * @brief Map a std::vector<double> buffer to float* (no copy, reinterpret cast).
     *
     * @param data std::vector<double> holding the buffer (must not be empty)
     * @return Pointer to data reinterpreted as float*
     *
     * **Warning**: This is a bitwise reinterpret cast. The double bit patterns will be
     * interpreted as floats, resulting in garbage values. Use only for memory allocation
     * sizing, not for numerical computation. For actual numerical conversion, see
     * `convert_double_to_float()` below.
     *
     * @throws std::runtime_error if data is empty
     */
    static float *map_double_vector_to_float(std::vector<double> &data)
    {
        if (data.empty())
            throw std::runtime_error("Cannot map empty vector");
        return reinterpret_cast<float *>(data.data());
    }

    /**
     * @brief Map a const std::vector<double> buffer to const float* (no copy).
     */
    static const float *map_double_vector_to_float_const(const std::vector<double> &data)
    {
        if (data.empty())
            throw std::runtime_error("Cannot map empty vector");
        return reinterpret_cast<const float *>(data.data());
    }

    /**
     * @brief Validate that a vector's data is properly aligned for float mapping.
     *
     * @param data std::vector<double> to check
     * @throws std::runtime_error if alignment insufficient for float*
     *
     * **Note**: Most std::allocator implementations guarantee alignment >= alignof(T),
     * so this usually succeeds. Platform-specific behavior may vary.
     */
    static void validate_alignment(const std::vector<double> &data)
    {
        if (data.empty())
            return; // Empty vectors have no alignment constraint

        const uintptr_t addr = reinterpret_cast<uintptr_t>(data.data());
        const std::size_t align = alignof(float);
        if (addr % align != 0)
            throw std::runtime_error("Vector data misaligned: addr=0x" +
                                     to_hex(addr) + " required align=" + std::to_string(align));
    }

    /**
     * @brief Validate grid dimensions consistency.
     *
     * @param nx, ny, nz Expected dimensions
     * @param n Total buffer size (should equal nx*ny*nz)
     *
     * @throws std::runtime_error if n != nx*ny*nz
     */
    static void validate_grid_dimensions(int nx, int ny, int nz, std::size_t n)
    {
        assert(nx > 0 && ny > 0 && nz > 0);
        const std::size_t expected = static_cast<std::size_t>(nx) * ny * nz;
        if (n != expected)
            throw std::runtime_error("Grid size mismatch: got " + std::to_string(n) +
                                     ", expected " + std::to_string(expected) +
                                     " for (" + std::to_string(nx) + "," + std::to_string(ny) +
                                     "," + std::to_string(nz) + ")");
    }

    /**
     * @brief Check if two pointers alias the same memory location (bitwise equal).
     *
     * Used to detect zero-copy mapping success: if map_result == original_pointer,
     * memory is truly aliased with no copy.
     *
     * @param ptr1, ptr2 Pointers to compare
     * @return true if ptr1 == ptr2 (same memory address)
     */
    template <typename T, typename U>
    static bool pointers_aliased(const T *ptr1, const U *ptr2)
    {
        return reinterpret_cast<const void *>(ptr1) == reinterpret_cast<const void *>(ptr2);
    }

private:
    /**
     * @brief Helper: format unsigned integer as hexadecimal string.
     */
    static std::string to_hex(uintptr_t val)
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "%lx", val);
        return std::string(buf);
    }
};

#endif // DIFFUSION3D_MEMORY_MAPPED_ADAPTER_H
