#pragma once

#include "diffusion3d_common.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

/**
 * @file vtk_export.h
 * @brief Export 3D scalar fields to VTK ImageData (.vti) for ParaView.
 */

#include <vector>

/**
 * @brief Export an arbitrary 3D field to VTK ImageData (.vti)
 * @param field Flat vector of size nx*ny*nz (z-fast order)
 * @param nx, ny, nz Grid dimensions
 * @param filename Output file path
 * @param varname Name of the variable (e.g., "delta")
 * @param viz_multiplier If not 1, writes physical `varname` plus `varname_viz` = field * multiplier
 *        (default active scalars: *_viz) so ParaView auto-range is not stuck at O(1) from early times.
 */
static void export_field_to_vti(
    const std::vector<double> &field,
    int nx, int ny, int nz,
    const std::string &filename,
    const std::string &varname = "delta",
    double viz_multiplier = 1.0)
{
    assert(viz_multiplier > 0.0);
    const bool dual = (viz_multiplier != 1.0);
    const std::string viz_name = varname + "_viz";

    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"1 1 1\">\n";
    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";
    f << "<PointData Scalars=\"" << (dual ? viz_name : varname) << "\">\n";
    f << "<DataArray type=\"Float64\" Name=\"" << varname << "\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                f << field[(z * ny + y) * nx + x] << " ";
            f << "\n";
        }
    f << "</DataArray>\n";
    if (dual)
    {
        f << "<DataArray type=\"Float64\" Name=\"" << viz_name << "\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
            {
                for (int x = 0; x < nx; ++x)
                    f << (field[(z * ny + y) * nx + x] * viz_multiplier) << " ";
                f << "\n";
            }
        f << "</DataArray>\n";
    }
    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
}

/**
 * @brief Export 3D diffusion field to VTK ImageData (.vti)
 *
 * This format can be opened directly in ParaView.
 */
static void export_to_vti(
    const Diffusion3DContext &ctx,
    const std::string &filename)
{
    std::ofstream f(filename);

    const int nx = ctx.nx;
    const int ny = ctx.ny;
    const int nz = ctx.nz;

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";

    std::ostringstream spacing;
    spacing << ctx.h << " " << ctx.h << " " << ctx.h;
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";

    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";

    f << "<PointData Scalars=\"u\">\n";
    f << "<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">\n";

    // VTK expects z-fast ordering consistent with ImageData
    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                f << ctx.u[ctx.idx(x, y, z)] << " ";
            }
            f << "\n";
        }
    }

    f << "</DataArray>\n";
    f << "</PointData>\n";

    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
}

/**
 * @brief Export CPU and GPU scalar fields on the same grid plus |cpu - gpu| for ParaView.
 *
 * One `.vti` with three point arrays: `cpu_u`, `gpu_u`, `abs_diff`.
 * In ParaView, use the field dropdown to color by each array.
 */
static void export_cpu_gpu_compare_to_vti(
    const std::vector<double> &cpu_field,
    const std::vector<double> &gpu_field,
    int nx, int ny, int nz,
    const std::string &filename)
{
    assert((int)cpu_field.size() == nx * ny * nz);
    assert((int)gpu_field.size() == nx * ny * nz);

    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"1 1 1\">\n";
    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";
    f << "<PointData>\n";

    f << "<DataArray type=\"Float64\" Name=\"cpu_u\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                f << cpu_field[(z * ny + y) * nx + x] << " ";
            f << "\n";
        }
    f << "</DataArray>\n";

    f << "<DataArray type=\"Float64\" Name=\"gpu_u\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                f << gpu_field[(z * ny + y) * nx + x] << " ";
            f << "\n";
        }
    f << "</DataArray>\n";

    f << "<DataArray type=\"Float64\" Name=\"abs_diff\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const size_t i = (size_t)((z * ny + y) * nx + x);
                const double d = cpu_field[i] - gpu_field[i];
                f << (d >= 0.0 ? d : -d) << " ";
            }
            f << "\n";
        }
    f << "</DataArray>\n";

    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
}