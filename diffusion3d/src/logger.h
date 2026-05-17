#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <vector>
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

/* ============================================================
 * STEP LOGGER (execution-level tracing)
 * ============================================================ */

struct StepLogger
{
    static void log(const std::string &stage)
    {
        std::cout << "[stage] " << stage << std::endl;
    }
};

#ifdef LOGGING
#define STEP_LOG(stage) StepLogger::log(stage)
#else
#define STEP_LOG(stage) ((void)0)
#endif

/* ============================================================
 * FIELD LOGGER (voxel-level, slice-level, CPU-only)
 * ============================================================ */

struct FieldLogger
{
    static void log_voxel(const std::string &label,
                          int x, int y, int z,
                          double value)
    {
        std::cout
            << "[field] " << label
            << "(" << x << "," << y << "," << z << ") = "
            << value << std::endl;
    }

    static void log_field_slice(const std::string &label,
                                int z,
                                const std::vector<double> &field,
                                int nx, int ny, int nz)
    {
        std::cout << "[field slice] " << label << " at z=" << z << std::endl;
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                double value = field[(z * ny + y) * nx + x];
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
    }
};

#ifdef LOGGING
#define FIELD_LOG(label, x, y, z, value) \
    FieldLogger::log_voxel(label, x, y, z, value)
#define FIELD_SLICE_LOG(label, z, field, nx, ny, nz) \
    FieldLogger::log_field_slice(label, z, field, nx, ny, nz)
#else
#define FIELD_LOG(label, x, y, z, value) ((void)0)
#define FIELD_SLICE_LOG(label, z, field, nx, ny, nz) ((void)0)
#endif

/* ============================================================
 * SLICE FIELD CAPTURE LOGGER (GPU-safe, time series volume slice)
 * ============================================================ */

#ifdef __CUDACC__
struct FieldSliceLogger
{
    int nx, ny, nz;
    int xMin, xMax, yMin, yMax, zMin, zMax;
    int steps;
    int sliceNx, sliceNy, sliceNz, sliceSize;
    std::vector<double> host_buffer;
    double *device_buffer = nullptr;

    FieldSliceLogger(int nx, int ny, int nz, int steps,
                     int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
        : nx(nx), ny(ny), nz(nz),
          xMin(xMin), xMax(xMax), yMin(yMin), yMax(yMax), zMin(zMin), zMax(zMax),
          steps(steps)
    {
        sliceNx = xMax - xMin + 1;
        sliceNy = yMax - yMin + 1;
        sliceNz = zMax - zMin + 1;
        sliceSize = sliceNx * sliceNy * sliceNz;
        host_buffer.resize(sliceSize * steps, 0.0);
        cudaMalloc(&device_buffer, sizeof(double) * sliceSize * steps);
        cudaMemset(device_buffer, 0, sizeof(double) * sliceSize * steps);
    }
    ~FieldSliceLogger()
    {
        if (device_buffer)
            cudaFree(device_buffer);
    }
    void fetch()
    {
        cudaMemcpy(host_buffer.data(), device_buffer, sizeof(double) * sliceSize * steps, cudaMemcpyDeviceToHost);
    }
    __host__ __device__ inline int idx(int t, int x, int y, int z) const
    {
        int sx = x - xMin;
        int sy = y - yMin;
        int sz = z - zMin;
        int sliceIndex = (sz * sliceNy + sy) * sliceNx + sx;
        return t * sliceSize + sliceIndex;
    }
    __host__ __device__ inline bool in_slice(int x, int y, int z) const
    {
        return (x >= xMin && x <= xMax &&
                y >= yMin && y <= yMax &&
                z >= zMin && z <= zMax);
    }
};
#else
struct FieldSliceLogger
{
    FieldSliceLogger(int, int, int, int, int, int, int, int, int, int) {}
    void fetch() {}
    int idx(int, int, int, int) const { return 0; }
    bool in_slice(int, int, int) const { return false; }
};
#endif

/* ============================================================
 * GPU SLICE WRITE MACRO
 * ============================================================ */

#ifdef LOGGING
/**
 * @brief Write a value to the logger's device buffer if the coordinate is within the slice region.
 * @param logger
 * @param t Timestep.
 * @param x
 * @param y
 * @param z
 * @param value
 */
#define SLICE_WRITE(logger, t, x, y, z, value) \
    if (logger.in_slice(x, y, z))              \
    {                                          \
        int i = logger.idx(t, x, y, z);        \
        logger.device_buffer[i] = (value);     \
    }
#else
#define SLICE_WRITE(logger, t, x, y, z, value) ((void)0)
#endif

#endif // LOGGER_H