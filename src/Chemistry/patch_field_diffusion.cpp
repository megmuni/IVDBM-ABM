#include "patch_field_diffusion.h"

#include "diffusion_algorithm_name.h"
#include "../Diffusion3D/core/make_multi_species_diffusion_engine.h"
#include "../Diffusion3D/core/multi_species_diffusion_settings.h"
#include "../Diffusion3D/core/multi_species_field_grid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

MultiSpeciesDiffusionSettings make_settings(const SpeciesRegistry &registry,
                                            DiffusionAlgorithm algorithm)
{
    MultiSpeciesDiffusionSettings settings;
    settings.algorithm = algorithm;
    for (SpeciesId id : registry.diffusing_species())
        settings.species_diffusivities[id] = registry.diffusivity(id);
    return settings;
}

} // namespace

struct PatchFieldDiffusion::Impl
{
    std::unique_ptr<MultiSpeciesFieldGrid> grid;
    std::unique_ptr<MultiSpeciesDiffusionEngine> engine;
    DiffusionAlgorithm engine_algorithm = DiffusionAlgorithm::ExplicitHeatEquation;
    std::vector<SpeciesId> species_ids;
    std::map<SpeciesId, std::vector<double>> initial_snapshot;

    void ensure_grid(int nx, int ny, int nz, double h, const SpeciesRegistry &registry)
    {
        const std::vector<SpeciesId> ids = registry.diffusing_species();
        const bool same_dims =
            grid && grid->nx_ == nx && grid->ny_ == ny && grid->nz_ == nz &&
            std::abs(grid->h_ - h) < 1e-12;
        const bool same_species =
            same_dims && species_ids.size() == ids.size() &&
            std::equal(species_ids.begin(), species_ids.end(), ids.begin());

        if (same_species)
            return;

        species_ids = ids;
        grid = std::make_unique<MultiSpeciesFieldGrid>(species_ids, nx, ny, nz, h);
        initial_snapshot.clear();
        const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;
        for (SpeciesId id : species_ids)
            initial_snapshot[id].assign(n, 0.0);
        engine.reset();
    }
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

void PatchFieldDiffusion::set_diffusion_algorithm(DiffusionAlgorithm algo)
{
    if (algorithm_ == algo)
        return;
    algorithm_ = algo;
    impl_->engine.reset();
}

DiffusionAlgorithm PatchFieldDiffusion::configured_algorithm() const
{
    return algorithm_;
}

const char *PatchFieldDiffusion::configured_algorithm_label() const
{
    return diffusion_algorithm_label(configured_algorithm());
}

const char *PatchFieldDiffusion::effective_algorithm_label() const
{
    if (impl_->engine)
        return diffusion_algorithm_label(impl_->engine->resolved_algorithm());
    return configured_algorithm_label();
}

void PatchFieldDiffusion::set_registry(const SpeciesRegistry &registry)
{
    if (registry.empty())
        throw std::invalid_argument("PatchFieldDiffusion: registry must not be empty");
    registry_ = &registry;
    impl_->engine.reset();
}

void PatchFieldDiffusion::diffuse_all_species(
    const std::map<SpeciesId, SpeciesDiffusionBuffers> &buffers,
    double tick_dt)
{
    if (registry_ == nullptr)
        throw std::runtime_error("PatchFieldDiffusion: set_registry() required before diffuse_all_species()");
    if (tick_dt <= 0.0)
        throw std::invalid_argument("PatchFieldDiffusion: tick_dt must be > 0");

    const std::vector<SpeciesId> ids = registry_->diffusing_species();
    if (ids.empty())
        throw std::runtime_error("PatchFieldDiffusion: registry has no diffusing species");

    for (SpeciesId id : ids)
    {
        const auto it = buffers.find(id);
        if (it == buffers.end())
            throw std::invalid_argument("PatchFieldDiffusion: missing buffers for species id");
        if (it->second.concentration == nullptr || it->second.diffused == nullptr)
            throw std::invalid_argument("PatchFieldDiffusion: null concentration or diffused buffer");
    }

    impl_->ensure_grid(nx_, ny_, nz_, h_, *registry_);

    const std::size_t n = static_cast<std::size_t>(nx_) * ny_ * nz_;

    for (SpeciesId id : ids)
    {
        const float *src = buffers.at(id).concentration;
        ScalarFieldGrid &field = impl_->grid->grid(id);
        std::vector<double> &snapshot = impl_->initial_snapshot[id];
        for (std::size_t i = 0; i < n; ++i)
        {
            snapshot[i] = static_cast<double>(src[i]);
            field.data_[i] = snapshot[i];
        }
    }

    const MultiSpeciesDiffusionSettings settings =
        make_settings(*registry_, algorithm_);
    const DiffusionAlgorithm resolved = ResolveDiffusionAlgorithm(settings);
    if (!impl_->engine || impl_->engine_algorithm != resolved)
    {
        impl_->engine = MakeMultiSpeciesDiffusionEngine(settings);
        impl_->engine_algorithm = resolved;
    }

    impl_->engine->configure_species_interval(*impl_->grid, settings, tick_dt);
    impl_->engine->advance_species_interval(*impl_->grid);

    for (SpeciesId id : ids)
    {
        const ScalarFieldGrid &field = impl_->grid->grid(id);
        float *dst = buffers.at(id).diffused;
        const std::vector<double> &snapshot = impl_->initial_snapshot.at(id);
        for (std::size_t i = 0; i < n; ++i)
            dst[i] = static_cast<float>(field.data_[i] - snapshot[i]);
    }
}
