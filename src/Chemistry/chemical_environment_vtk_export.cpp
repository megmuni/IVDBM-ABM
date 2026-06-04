#include "chemical_environment_vtk_export.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace
{

void write_vti_array(std::ofstream &f, const std::string &name,
                     const float *grid, int nx, int ny, int nz,
                     double viz_multiplier)
{
    f << "<DataArray type=\"Float64\" Name=\"" << name
      << "\" format=\"ascii\">\n";
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int idx = x + nx * (y + ny * z);
                f << static_cast<double>(grid[idx]) << " ";
            }
            f << "\n";
        }
    f << "</DataArray>\n";

    if (viz_multiplier != 1.0)
    {
        const std::string viz_name = name + "_viz";
        f << "<DataArray type=\"Float64\" Name=\"" << viz_name
          << "\" format=\"ascii\">\n";
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
            {
                for (int x = 0; x < nx; ++x)
                {
                    const int idx = x + nx * (y + ny * z);
                    f << static_cast<double>(grid[idx]) * viz_multiplier
                      << " ";
                }
                f << "\n";
            }
        f << "</DataArray>\n";
    }
}

} // namespace

std::vector<ChemicalEnvironmentVtkChannel>
default_vtk_channels(const ChemicalEnvironment &env,
                     const ChemicalEnvironmentVtkExportOptions &options)
{
    std::vector<ChemicalEnvironmentVtkChannel> channels;

    if (options.export_concentrations || options.export_delta_channels)
    {
        for (SpeciesId id : env.registry().diffusing_species())
        {
            const SpeciesDescriptor &desc = env.registry().descriptor(id);
            if (options.export_concentrations)
            {
                channels.push_back(
                    {desc.name, desc.concentration_channel});
            }
            if (options.export_delta_channels)
            {
                channels.push_back(
                    {"d" + desc.name, desc.diffused_channel});
            }
        }
    }

    if (options.export_chemotaxis && env.chemotaxis_channel_index() >= 0)
    {
        channels.push_back({options.chemotaxis_array_name,
                            env.chemotaxis_channel_index()});
    }

    for (const ChemicalEnvironmentVtkChannel &extra : options.extra_channels)
        channels.push_back(extra);

    return channels;
}

bool export_chemical_environment_to_vti(
    const ChemicalEnvironment &env, const std::string &filename,
    const ChemicalEnvironmentVtkExportOptions &options)
{
    if (env.channel_count() <= 0)
        return false;

    const int nx = env.grid_nx();
    const int ny = env.grid_ny();
    const int nz = env.grid_nz();
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return false;

    double h = options.grid_spacing_mm;
    if (h <= 0.0)
        h = env.grid_spacing_mm();
    if (h <= 0.0)
        throw std::invalid_argument(
            "export_chemical_environment_to_vti: grid_spacing_mm must be > 0");

    const std::vector<ChemicalEnvironmentVtkChannel> arrays =
        default_vtk_channels(env, options);
    if (arrays.empty())
        return false;

    for (const ChemicalEnvironmentVtkChannel &ch : arrays)
    {
        if (ch.channel_index < 0 || ch.channel_index >= env.channel_count())
            return false;
        if (ch.array_name.empty())
            return false;
    }

    std::ostringstream spacing;
    spacing << h << " " << h << " " << h;

    std::ofstream f(filename);
    if (!f)
        return false;

    const bool dual = (options.viz_multiplier != 1.0);
    const std::string active =
        dual ? arrays.front().array_name + "_viz" : arrays.front().array_name;

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "<ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 "
      << nz - 1 << "\" "
      << "Origin=\"0 0 0\" "
      << "Spacing=\"" << spacing.str() << "\">\n";
    f << "<Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1
      << "\">\n";
    f << "<PointData Scalars=\"" << active << "\">\n";

    for (const ChemicalEnvironmentVtkChannel &ch : arrays)
    {
        const float *grid = env.channel_grid(ch.channel_index);
        if (grid == nullptr)
            return false;
        write_vti_array(f, ch.array_name, grid, nx, ny, nz,
                        options.viz_multiplier);
    }

    f << "</PointData>\n";
    f << "</Piece>\n";
    f << "</ImageData>\n";
    f << "</VTKFile>\n";
    return true;
}

std::string format_chemical_environment_vti_path(const std::string &directory,
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
