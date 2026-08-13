# validate_rcwa.jl — confirm the RCWA engine before using it for the metasurface workflow.
#
# Three independent checks:
#   (1) ENERGY CONSERVATION: a lossless dielectric pattern must satisfy R + T = 1.
#   (2) RCWA ↔ TMM cross-check: a uniform (full-coverage) VO2|Al stack must match the
#       transfer-matrix result to near machine precision.
#   (3) CONVERGENCE: the spectrum must converge as the RCWA truncation order N increases.
#
# Run:  julia --project=. -t auto metasurface_ir/validate_rcwa.jl

using GESTALT
using Plots
using LinearAlgebra
include(joinpath(@__DIR__, "metasurface_ir.jl"))
BLAS.set_num_threads(1)

const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)
fig(name) = joinpath(FIGDIR, name)

p = 5e-6
air = HomogeneousLayer(NKMaterial(1.0), Inf)
vo2i = vo2_insulating_material()
al = aluminum_material()

# =====================================================================================
println("(1) Energy conservation on a lossless silicon disk pattern (R + T = 1) ...")
si = NKMaterial(3.46)
disk = circle_image(0.5, 0.5, 0.32)
ind = Float64.(disk .> 0)
pat = fourier_pattern_from_indicator(ind, p, si, NKMaterial(1.0); mmax=16, nmax=16)
layer = PeriodicLayer(si, 2e-6, pat)
grid8 = RCWAGrid((2π / p, 0.0), (0.0, 2π / p); Nx=8, Ny=8)
lossless = [air, layer, air]
λs1 = range(8e-6, 14e-6; length=16)
_, _, A1 = sweep_spectrum(lossless, grid8, λs1; unpolarized=false)
println("    max |R+T−1| = max |A| = $(maximum(abs.(A1)))   (expect ≈ 1e-12)")

# =====================================================================================
println("\n(2) RCWA ↔ TMM cross-check for uniform VO2(4 µm)|Al(200 nm) ...")
full = fill(UInt8(0xff), 64, 64)
grid6, rcwa_layers = build_stack(full, vo2i, al; N=6)
vo2l = HomogeneousLayer(vo2i, 4e-6)
all  = HomogeneousLayer(al, 200e-9)
tmm_layers = [air, vo2l, all, air]
for λ in (8e-6, 10e-6, 12e-6, 14e-6)
    f = c₀ / λ
    rr = solve(rcwa_layers, grid6, f; θ=0.0, pol_TE=1.0)
    tt = solve(tmm_layers, f; θ=0.0)
    println("    λ=$(λ*1e6) µm  RCWA R=$(round(rr.Rtot, digits=6)) T=$(round(rr.Ttot, digits=9))" *
            "  |  TMM R=$(round(tt.R, digits=6)) T=$(round(tt.T, digits=9))" *
            "  ΔR=$(abs(rr.Rtot - tt.R))  ΔT=$(abs(rr.Ttot - tt.T))")
end

# =====================================================================================
println("\n(3) Convergence of absorptance with RCWA truncation N ...")
shape = rect_image(0.2, 0.2, 0.5, 0.35)
λs3 = range(8e-6, 14e-6; length=24)
i10 = argmin(abs.(λs3 .- 10e-6))
i13 = argmin(abs.(λs3 .- 13e-6))
pN = plot(; xlabel="λ (µm)", ylabel="absorptance A", title="Convergence vs truncation N",
          legend=:topright)
for N in (2, 4, 6, 8)
    g, L = build_stack(shape, vo2i, al; N=N)
    _, _, A = sweep_spectrum(L, g, λs3; unpolarized=false)
    println("    N=$N  A(10µm)=$(round(A[i10], digits=4))  A(13µm)=$(round(A[i13], digits=4))")
    plot!(pN, λs3 .* 1e6, A; lw=2, label="N=$N")
end
save_plot(pN, fig("validate_convergence.png"))

println("\nDone.")
