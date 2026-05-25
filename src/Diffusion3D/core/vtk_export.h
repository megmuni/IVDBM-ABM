#ifndef DIFFUSION3D_VTK_EXPORT_H
#define DIFFUSION3D_VTK_EXPORT_H

/**
 * @file vtk_export.h
 * @brief Export 3D scalar fields to VTK ImageData (.vti) for ParaView.
 *
 * Adapted from diffusion3d/src/vtk_export.h. Uses ScalarFieldGrid / MultiSpeciesFieldGrid
 * (row-major layout: idx = x + nx*(y + ny*z), compatible with VTK ImageData point order).
 */

#include "multi_species_field_grid.h"
#include "scalar_field_grid.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @brief Export an arbitrary 3D field to VTK ImageData (.vti).
 *
 * @param field Flat vector of size nx*ny*nz (same layout as ScalarFieldGrid::data_)
 * @param nx, ny, nz Grid dimensions
 * @param h Uniform grid spacing written to VTK Spacing (same units as simulation)
 * @param filename Output file path
 * @param varname Name of the scalar array (e.g. "concentration")
 * @param viz_multiplier If not 1, also writes `varname_viz` = field * multiplier for ParaView auto-range
 */
inline void export_field_to_vti(
    const std::vector<double> &field,
    int nx, int ny, int nz,
    double h,
    const std::string &filename,
    const std::string &varname = "u",
    double viz_multiplier = 1.0)
{
    assert(nx > 0 && ny > 0 && nz > 0);
    assert(h > 0.0);
    assert(static_cast<int>(field.size()) == nx * ny * nz);
    assert(viz_multiplier > 0.0);

    const bool dual = (viz_multiplier != 1.0);
    const std::string viz_name = varname + "_viz";

    std::ostringstream spacing;
    spacing << h << " " << h << " " << h;

    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";
    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";
    f << "<PointData Scalars=\"" << (dual ? viz_name : varname) << "\">\n";
    f << "<DataArray type=\"Float64\" Name=\"" << varname << "\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                f << field[x + nx * (y + ny * z)] << " ";
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
                    f << (field[x + nx * (y + ny * z)] * viz_multiplier) << " ";
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
 * @brief Export one species grid to `.vti` (opens directly in ParaView).
 */
inline void export_scalar_field_to_vti(
    const ScalarFieldGrid &grid,
    double h,
    const std::string &filename,
    const std::string &varname = "u",
    double viz_multiplier = 1.0)
{
    export_field_to_vti(grid.data_, grid.nx_, grid.ny_, grid.nz_, h, filename, varname, viz_multiplier);
}

/**
 * @brief Export all species in one `.vti` with one point array per species.
 *
 * Array names: `species_<id>` (e.g. `species_0`, `species_1`). Active scalars: first species.
 */
inline void export_multi_species_field_to_vti(
    const MultiSpeciesFieldGrid &grid,
    const std::string &filename)
{
    assert(!grid.species().empty());

    const int nx = grid.nx_;
    const int ny = grid.ny_;
    const int nz = grid.nz_;
    const double h = grid.h_;

    std::ostringstream spacing;
    spacing << h << " " << h << " " << h;

    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";
    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";

    const auto &species = grid.species();
    f << "<PointData Scalars=\"species_" << species.front() << "\">\n";

    for (SpeciesId id : species)
    {
        const ScalarFieldGrid &g = grid.grid(id);
        const std::string name = "species_" + std::to_string(id);
        f << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
            {
                for (int x = 0; x < nx; ++x)
                    f << g.data_[g.idx(x, y, z)] << " ";
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
 * @brief Export CPU and GPU fields plus |cpu - gpu| in one `.vti` for backend comparison in ParaView.
 */
inline void export_cpu_gpu_compare_to_vti(
    const std::vector<double> &cpu_field,
    const std::vector<double> &gpu_field,
    int nx, int ny, int nz,
    double h,
    const std::string &filename)
{
    assert(static_cast<int>(cpu_field.size()) == nx * ny * nz);
    assert(static_cast<int>(gpu_field.size()) == nx * ny * nz);

    std::ostringstream spacing;
    spacing << h << " " << h << " " << h;

    std::ofstream f(filename);
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";
    f << "<Piece Extent=\"0 " << nx - 1
      << " 0 " << ny - 1
      << " 0 " << nz - 1 << "\">\n";
    f << "<PointData Scalars=\"cpu_u\">\n";

    auto write_array = [&](const char *name, auto value_at) {
        f << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
            {
                for (int x = 0; x < nx; ++x)
                    f << value_at(x, y, z) << " ";
                f << "\n";
            }
        f << "</DataArray>\n";
    };

    const auto lin = [&](int x, int y, int z) { return x + nx * (y + ny * z); };

    write_array("cpu_u", [&](int x, int y, int z) { return cpu_field[lin(x, y, z)]; });
    write_array("gpu_u", [&](int x, int y, int z) { return gpu_field[lin(x, y, z)]; });
    write_array("abs_diff", [&](int x, int y, int z) {
        const double d = cpu_field[lin(x, y, z)] - gpu_field[lin(x, y, z)];
        return d >= 0.0 ? d : -d;
    });

    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
}

#endif // DIFFUSION3D_VTK_EXPORT_H
