#include "world_vtk_export.h"

#include "Usr_World/biomaterialWorld.h"
#include "World.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace
{

void write_spacing(std::ostringstream &spacing, double sx, double sy, double sz)
{
    spacing << sx << " " << sy << " " << sz;
}

bool resolve_spacing(const World &world, const WorldVtkExportOptions &options,
                     double &sx, double &sy, double &sz)
{
    if (options.grid_spacing_mm > 0.0)
    {
        sx = sy = sz = options.grid_spacing_mm;
        return true;
    }

    sx = world.dx;
    sy = world.dy;
    sz = world.dz;
    return sx > 0.0 && sy > 0.0 && sz > 0.0;
}

bool ensure_directory(const std::string &directory)
{
    if (directory.empty())
        return false;

    std::error_code ec;
    std::filesystem::create_directories(directory, ec);
    return !ec;
}

} // namespace

std::string format_world_vti_path(const std::string &directory,
                                  const std::string &basename,
                                  int step_index,
                                  int zero_pad)
{
    std::ostringstream path;
    path << directory;
    if (!directory.empty() && directory.back() != '/')
        path << "/";
    path << basename << "_t" << std::setw(zero_pad) << std::setfill('0')
         << step_index << ".vti";
    return path.str();
}

bool export_world_patches_to_vti(const World &world, const std::string &filename,
                                 const WorldVtkExportOptions &options)
{
    const int nx = world.nx;
    const int ny = world.ny;
    const int nz = world.nz;
    if (nx <= 0 || ny <= 0 || nz <= 0 || world.worldPatch == nullptr)
        return false;

    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    if (!resolve_spacing(world, options, sx, sy, sz))
        return false;

    std::ostringstream spacing;
    write_spacing(spacing, sx, sy, sz);

    std::ofstream f(filename);
    if (!f)
        return false;

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 "
      << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";
    f << "<Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1
      << "\">\n";
    f << "<PointData Scalars=\"occupied\">\n";
    f << "<DataArray type=\"Int32\" Name=\"occupied\" format=\"ascii\">\n";

    for (int z = 0; z < nz; ++z)
    {
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int idx = x + nx * (y + ny * z);
                f << (world.worldPatch[idx].isOccupied() ? 1 : 0) << " ";
            }
            f << "\n";
        }
    }

    f << "</DataArray>\n";
    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
    return true;
}

bool export_world_timestep(const BMWorld &world, int step_index,
                           const std::string &directory,
                           const WorldVtkExportOptions &options)
{
    if (!ensure_directory(directory))
        return false;

    const std::string patches_path =
        format_world_vti_path(directory, "patches", step_index);
    if (!export_world_patches_to_vti(world, patches_path, options))
        return false;

    const ChemicalEnvironment *chem_env = world.chemical_environment();
    if (chem_env == nullptr)
        return true;

    ChemicalEnvironmentVtkExportOptions chem_options = options.chem_options;
    if (chem_options.grid_spacing_mm <= 0.0 && options.grid_spacing_mm > 0.0)
        chem_options.grid_spacing_mm = options.grid_spacing_mm;

    const std::string chem_path =
        format_world_vti_path(directory, "chem", step_index);
    export_chemical_environment_to_vti(*chem_env, chem_path, chem_options);
    return true;
}
