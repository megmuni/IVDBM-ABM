#include "PatchFieldDiffusion.h"

#include <stdexcept>

struct PatchFieldDiffusion::Impl
{
};

PatchFieldDiffusion::PatchFieldDiffusion(int nx, int ny, int nz, double grid_spacing)
    : impl_(new Impl()), nx_(nx), ny_(ny), nz_(nz), h_(grid_spacing)
{
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0)
        throw std::invalid_argument("PatchFieldDiffusion: grid dimensions must be positive");
    if (h_ <= 0.0)
        throw std::invalid_argument("PatchFieldDiffusion: grid spacing must be > 0");
}

PatchFieldDiffusion::~PatchFieldDiffusion() = default;

void PatchFieldDiffusion::set_registry(const SpeciesRegistry &registry)
{
    if (registry.empty())
        throw std::invalid_argument("PatchFieldDiffusion: registry must not be empty");
    registry_ = &registry;
}

void PatchFieldDiffusion::diffuse_all_species(
    const std::map<SpeciesId, SpeciesDiffusionBuffers> &buffers,
    double tick_dt)
{
    if (registry_ == nullptr)
        throw std::runtime_error("PatchFieldDiffusion: set_registry() required before diffuse_all_species()");
    if (tick_dt <= 0.0)
        throw std::invalid_argument("PatchFieldDiffusion: tick_dt must be > 0");

    const std::size_t n = static_cast<std::size_t>(nx_) * ny_ * nz_;
    for (SpeciesId id : registry_->diffusing_species())
    {
        const auto it = buffers.find(id);
        if (it == buffers.end())
            throw std::invalid_argument("PatchFieldDiffusion: missing buffers for species id");
        if (it->second.concentration == nullptr || it->second.diffused == nullptr)
            throw std::invalid_argument("PatchFieldDiffusion: null concentration or diffused buffer");
        (void)n;
        (void)it;
    }

    throw std::runtime_error("PatchFieldDiffusion::diffuse_all_species not implemented (Phase II.3)");
}
