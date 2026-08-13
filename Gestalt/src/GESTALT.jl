"""
    GESTALT

A General Electromagnetic Solver for Arbitrarily Twisted and Layered photonic Topologies.

GESTALT models electromagnetic wave propagation through generalised stacked photonic
structures. It combines two complementary engines:

  * **Rigorous Coupled-Wave Analysis (RCWA)** — a full-wave, Fourier-space method that
    returns exact scattering matrices for periodic (and homogeneous) layers under
    plane-wave illumination, including polarisation and diffraction orders.
  * **Transfer-Matrix Method (TMM)** — the classic 2×2 stratified-media method for
    laterally homogeneous multilayer stacks, used as a fast baseline and a reference.

The building blocks are organised as:

  * `Material.jl` — dispersive bulk material models (`AbstractMaterial` and `𝜀̃`, `ñ`).
  * `Layer.jl`    — layers built on a material + thickness: homogeneous, periodic
                    (RCWA-compatible) and aperiodic in-plane patterns.
  * `TMM.jl`      — transfer-matrix solver for homogeneous stacks.
  * `RCWA.jl`     — rigorous coupled-wave analysis solver and scattering matrices.

Both engines expose a common `solve` entry point (dispatched on the argument list) and a
consistent SI unit convention (frequency `f` in Hz, lengths in metres, time dependence
e^{-iωt}).

# Quick start
```julia
using GESTALT

# Homogeneous multilayer (TMM): air | quarter-wave film | glass
air   = HomogeneousLayer(NKMaterial(1.0), Inf)
film  = HomogeneousLayer(NKMaterial(1.22), 112e-9)
glass = HomogeneousLayer(NKMaterial(1.5), Inf)
res   = solve([air, film, glass], c₀ / 550e-9)      # res.R, res.T, res.A

# 2D photonic-crystal slab (RCWA)
grid  = RCWAGrid((2π/1e-6, 0.0), (0.0, 2π/1e-6); Nx=4, Ny=4)
slab  = PeriodicLayer(NKMaterial(3.46), 0.5e-6, pattern)   # pattern :: FourierPeriodicPattern
out   = solve([air, slab, air], grid, f)            # out.Rtot, out.Ttot, out.R, out.T
```
"""
module GESTALT

# Implementation files. Each defines its types/functions in this module's scope; the
# include-guards inside them are no-ops here because everything is loaded in order.
include("Material.jl")
include("Layer.jl")
include("TMM.jl")
include("RCWA.jl")
include("FDFD.jl")
include("RCWA4D.jl")
include("CMT.jl")
include("Visualize.jl")

# ---- Physical constants ----
export c₀

# ---- Materials (Material.jl) ----
export AbstractMaterial,
       NKMaterial, FunctionalNKMaterial, FunctionalEpsilonMaterial,
       LorentzMaterial, DrudeMaterial, DebyeMaterial,
       SampledDataMaterial,
       𝜀̃, ñ

# ---- Layers and in-plane patterns (Layer.jl) ----
export AbstractLayer,
       AbstractHomogeneousLayer, HomogeneousLayer,
       AbstractPeriodicLayer, PeriodicLayer,
       AbstractAperiodicLayer, AperiodicLayer,
       AbstractPeriodicPattern, FourierPeriodicPattern, FunctionalPeriodicPattern,
       PixelPeriodicPattern,
       AbstractAperiodicPattern, FunctionalAperiodicPattern, GridAperiodicPattern,
       material, thickness, pattern,
       ε̂, g_vectors, ε̃_inplane,
       𝜀̃_material, ñ_material

# ---- Transfer-matrix method (TMM.jl) ----
export AbstractPolarization, TE, TM,
       TMMResult,
       reflectance, transmittance, absorptance

# ---- Rigorous coupled-wave analysis (RCWA.jl) ----
export RCWAGrid, nharmonics, convolution_matrix,
       SMatrix, star, smatrix_k, resonance_indicator,
       RCWAResult

# ---- Twisted multilayer RCWA4D (RCWA4D.jl) ----
export TwistedGrid, ndevicelayers, joint_convolution_matrix,
       smatrix_k_twisted, solve_twisted

# ---- Tight-binding coupled-mode theory (CMT.jl) ----
export cmt_matrix, eigenfrequencies, resonances, linewidths, qfactors,
       CMTCoupling, bloch_matrix, bands,
       extract_resonance, extract_coupling

# ---- Common solver (a TMM method and an RCWA method of the same generic) ----
export solve

# ---- Visualization (Visualize.jl) ----
export save_plot, plot_material, plot_rta,
       plot_reciprocal_basis, plot_pattern,
       plot_diffraction_orders, plot_bands,
       plot_stack, plot_band_map

end # module GESTALT
