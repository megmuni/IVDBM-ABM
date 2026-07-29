#ifndef IVDBM_WORLD_VTK_EXPORT_H
#define IVDBM_WORLD_VTK_EXPORT_H

/**
 * @file world_vtk_export.h
 * @brief Export @ref World patch occupancy and @ref BMWorld chemistry to VTK
 *        ImageData (.vti) for ParaView.
 */

#include "../Chemistry/chemical_environment_vtk_export.h"

#include <string>

class World;

/** Controls VTK metadata for World export. */
struct WorldVtkExportOptions
{
    /**
     * Uniform patch spacing in mm for VTK Spacing when > 0.
     * If <= 0, uses @c world.dx, @c world.dy, @c world.dz.
     */
    double grid_spacing_mm = 0.0;

    /** Passed through to @ref export_chemical_environment_to_vti. */
    ChemicalEnvironmentVtkExportOptions chem_options{};
};

/**
 * @brief Write one `.vti` with a single @c occupied point array (0/1).
 *
 * Works for any @ref World subclass; reads @c worldPatch via @ref Patch::isOccupied().
 */
bool export_world_patches_to_vti(const World &world, const std::string &filename,
                                 const WorldVtkExportOptions &options = {});

class BMWorld;

/**
 * @brief Write @c patches_t<step>.vti and, when allocated, @c chem_t<step>.vti
 *        under @p directory.
 */
bool export_world_timestep(const BMWorld &world, int step_index,
                           const std::string &directory,
                           const WorldVtkExportOptions &options = {});

/** Path helper for time-series export, e.g. @c dir/patches_t00003.vti . */
std::string format_world_vti_path(const std::string &directory,
                                  const std::string &basename,
                                  int step_index,
                                  int zero_pad = 5);

#endif
