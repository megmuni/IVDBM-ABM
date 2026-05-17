# Diffusion3D Module

This module implements multi-species scalar field diffusion with zero-copy memory mapping and shared pointer support. It is designed for integration with ABM frameworks and supports extensible PDE stepping strategies.

## Structure

- `core/` - Core diffusion engine, field grids, settings, and stepper implementations
- `adapters/` - Memory mapping and type adapters
- `tests/` - Unit and integration tests

## Key Types

- `MultiSpeciesFieldGrid`
- `MultiSpeciesDiffusionSettings`
- `MultiSpeciesDiffusionEngine`
- `ExplicitMultiSpeciesHeatStepper`
- `MakeMultiSpeciesDiffusionEngine()`
- `MemoryMappedAdapter`

## Build

Add this module to your CMake project and enable with `BUILD_DIFFUSION3D`.
