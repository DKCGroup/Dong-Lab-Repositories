# run_shapes.jl — IR metasurface (8–14 µm) RCWA simulation for randomly-generated top patterns.
#
# Structure (period p = 5 µm):
#     incidence: air
#     top:       VO2 pattern, height from the grayscale image (gray/255 × 4 µm)
#     bottom:    200 nm aluminium mirror (opaque ⇒ T ≈ 0, so A = 1 − R)
# VO2 is modelled in both its insulating (Lorentz) and metallic (Drude) phases; measured n,k
# (SampledDataMaterial) can be swapped in later.
#
# Run:  julia --project=. -t auto metasurface_ir/run_shapes.jl

using GESTALT
using Plots
using LinearAlgebra
using Random
include(joinpath(@__DIR__, "metasurface_ir.jl"))
BLAS.set_num_threads(1)

const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)
fig(name) = joinpath(FIGDIR, name)
println("Writing figures to: ", FIGDIR, "   (threads = ", Threads.nthreads(), ")")

# ---- Parameters -------------------------------------------------------------------
p   = 5e-6                       # period (m)
max_t = 4e-6                     # grayscale 255 → 4 µm
λs  = range(8e-6, 14e-6; length=256)   # ≥ 200 frequency samples
N   = 6                          # RCWA truncation (NH = (2N+1)² = 169); fine Fourier grid.
                                 #   validate_rcwa.jl shows N ≥ 6 is converged to <0.5%.
                                 #   N = 8 (NH = 289) is available for publication runs.

vo2i = vo2_insulating_material()
vo2m = vo2_metallic_material()
al   = aluminum_material()

# ---- Random representative shapes (seeded for reproducibility) --------------------
Random.seed!(2024)
shapes = [
    ("rectangle", rect_image(0.15, 0.25, 0.60, 0.40)),
    ("circle",    circle_image(0.50, 0.50, 0.34)),
    ("triangle",  triangle_image((0.12, 0.12), (0.88, 0.18), (0.42, 0.82))),
    ("cross",     cross_image(0.20, 0.20, 0.60, 0.60, 0.22)),
]

# ---- Per-shape simulation + figures ----------------------------------------------
function thickness_map_plot(img)
    h = thickness_profile(img, max_t) .* 1e6                       # µm
    Nimg = size(img, 1)
    xs = range(0, p; length=Nimg) .* 1e6
    ys = range(0, p; length=Nimg) .* 1e6
    heatmap(xs, ys, permutedims(h); aspect_ratio=:equal, xlabel="x (µm)", ylabel="y (µm)",
            title="pattern thickness (µm)", color=:viridis, colorbar_title="µm",
            clims=(0, max_t * 1e6))
end

function spectrum_plot(Ri, Ai, Rm, Am, Ti, Tm)
    λµm = λs .* 1e6
    pl = plot(λµm, Ri; lw=2, label="R  ins", color=:steelblue,
              xlabel="λ (µm)", ylabel="efficiency", ylims=(0, 1),
              title="8–14 µm R / T / A", legend=:topright)
    plot!(pl, λµm, Ai; lw=2, ls=:dash, label="A  ins", color=:steelblue)
    plot!(pl, λµm, Rm; lw=2, label="R  met", color=:indianred)
    plot!(pl, λµm, Am; lw=2, ls=:dash, label="A  met", color=:indianred)
    plot!(pl, λµm, Ti; lw=1.5, ls=:dot, label="T  ins", color=:seagreen)
    plot!(pl, λµm, Tm; lw=1.5, ls=:dot, label="T  met", color=:seagreen)
    pl
end

summary = NamedTuple[]
for (name, img) in shapes
    println("\nshape: $(name)   (fill $(round(100*sum(img .> 0)/length(img), digits=1))%)")
    grid, L_ins = build_stack(img, vo2i, al; p=p, max_t=max_t, N=N)
    _,    L_met = build_stack(img, vo2m, al; p=p, max_t=max_t, N=N)

    t = @elapsed begin
        Ri, Ti, Ai = sweep_spectrum(L_ins, grid, λs)
        Rm, Tm, Am = sweep_spectrum(L_met, grid, λs)
    end
    push!(summary, (name=name, Ri=Ri, Ai=Ai, Rm=Rm, Am=Am))
    println("    swept $(length(λs)) λ × 2 phases (TE+TM) in $(round(t, digits=1)) s")

    # sanity: energy conservation by construction (A = 1 − R − T), and physical bounds
    println("    R: $(round(minimum(Ri), digits=3))–$(round(maximum(Ri), digits=3)) (ins) / " *
            "$(round(minimum(Rm), digits=3))–$(round(maximum(Rm), digits=3)) (met)")
    println("    A: $(round(minimum(Ai), digits=3))–$(round(maximum(Ai), digits=3)) (ins) / " *
            "$(round(minimum(Am), digits=3))–$(round(maximum(Am), digits=3)) (met)")

    save_plot(plot(thickness_map_plot(img), spectrum_plot(Ri, Ai, Rm, Am, Ti, Tm);
                   layout=(1, 2), size=(1150, 450), left_margin=5Plots.mm,
                   bottom_margin=5Plots.mm, plot_title=name),
              fig("shape_$(name).png"))
end

# ---- Summary: absorption of all shapes -------------------------------------------
psum = plot(; xlabel="λ (µm)", ylabel="absorptance A", ylims=(0, 1),
            title="Absorption of all shapes (solid = ins, dashed = met)", legend=:topright)
palette = [:steelblue, :indianred, :seagreen, :darkorange]
for (k, s) in enumerate(summary)
    plot!(psum, λs .* 1e6, s.Ai; lw=2, ls=:solid, label="$(s.name) ins", color=palette[k])
    plot!(psum, λs .* 1e6, s.Am; lw=2, ls=:dash,  label="$(s.name) met", color=palette[k])
end
save_plot(psum, fig("summary_absorption.png"))

println("\nDone. Wrote per-shape figures + summary to ", FIGDIR)
