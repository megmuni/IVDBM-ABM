#ifndef IVDBM_CHEMICAL_ENVIRONMENT_VTK_EXPORT_H
#define IVDBM_CHEMICAL_ENVIRONMENT_VTK_EXPORT_H

/**
 * @file chemical_environment_vtk_export.h
 * @brief Export @ref ChemicalEnvironment channel grids to VTK ImageData (.vti) for ParaView.
 */

#include "chemical_environment.h"

#include <string>
#include <vector>

/** One point array to write into a multi-channel .vti file. */
struct ChemicalEnvironmentVtkChannel
{
    std::string array_name;
    int channel_index = -1;
};

/** Controls which grids are written and VTK metadata. */
struct ChemicalEnvironmentVtkExportOptions
{
    /**
     * Patch spacing in mm for VTK Spacing. If <= 0, uses
     * @ref ChemicalEnvironment::grid_spacing_mm() (must be > 0).
     */
    double grid_spacing_mm = 0.0;

    /**
     * If not 1, also writes @c array_name_viz = field * multiplier for each array
     * (helps ParaView auto-range on small concentrations).
     */
    double viz_multiplier = 1.0;

    /** Export each registered species concentration channel (default on). */
    bool export_concentrations = true;

    /** Export each species diffused/delta channel (e.g. dTNF). */
    bool export_delta_channels = false;

    /** Export the chemotaxis channel when allocated. */
    bool export_chemotaxis = false;

    /** Name for chemotaxis array when @c export_chemotaxis is true. */
    std::string chemotaxis_array_name = "chemotaxis";

    /** Additional channel grids (appended after built-in selections). */
    std::vector<ChemicalEnvironmentVtkChannel> extra_channels;
};

/**
 * @brief Build the list of VTK point arrays from @p options and @p env registry/config.
 */
std::vector<ChemicalEnvironmentVtkChannel>
default_vtk_channels(const ChemicalEnvironment &env,
                     const ChemicalEnvironmentVtkExportOptions &options);

/**
 * @brief Write one VTK ImageData file with all selected channel arrays.
 *
 * Grid layout matches the patch index order: @c idx = x + nx * (y + ny * z).
 *
 * @return false if the file could not be written or channels are not allocated.
 */
bool export_chemical_environment_to_vti(
    const ChemicalEnvironment &env, const std::string &filename,
    const ChemicalEnvironmentVtkExportOptions &options = {});

/**
 * @brief Path helper for time-series export, e.g. @c dir/chem_t00003.vti .
 */
std::string format_chemical_environment_vti_path(const std::string &directory,
                                                 const std::string &basename,
                                                 int step_index,
                                                 int zero_pad = 5);

#endif
