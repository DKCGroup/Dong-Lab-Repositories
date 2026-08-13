using Plots

# Visualize — lightweight plotting helpers for GESTALT results, built on Plots.jl. Every
# function returns the Plots plot object so it can be further customised, displayed, or written
# to disk with `save_plot`. The helpers cover the quantities the framework already produces:
# material dispersion, R/T/A spectra, RCWA diffraction orders, the reciprocal Fourier module
# (single-lattice or twisted joint basis), real-space permittivity patterns, and CMT bands.
#
# Loaded last in the package so all result/grid/material types are already defined; the guards
# below let this file also be `include`d standalone.
if !isdefined(@__MODULE__, :TMMResult)
    include(joinpath(@__DIR__, "TMM.jl"))
end
if !isdefined(@__MODULE__, :TwistedGrid)
    include(joinpath(@__DIR__, "RCWA4D.jl"))
end
if !isdefined(@__MODULE__, :bands)
    include(joinpath(@__DIR__, "CMT.jl"))
end


# ---- Small helpers ----
# Take the requested scalar component (:real, :imag, or :abs) of a complex value.
_partval(v, part::Symbol) =
    part === :real ? real(v) :
    part === :imag ? imag(v) :
    part === :abs  ? abs(v)  :
    error("part must be :real, :imag, or :abs, got $part")

# Reflectance / transmittance / absorptance triple from either result type.
_rta(r::TMMResult)  = (r.R, r.T, r.A)
_rta(r::RCWAResult) = (r.Rtot, r.Ttot, r.A)

# Write a plot to disk (PNG/PDF/SVG inferred from the extension) and return the path.
save_plot(plt, path::AbstractString) = (Plots.savefig(plt, path); path)


# ---- Material dispersion ----
# Plot a material's complex permittivity (`quantity=:epsilon`) or refractive index
# (`quantity=:index`) versus frequency (THz). `fs` are frequencies in Hz.
function plot_material(material::AbstractMaterial, fs::AbstractVector{<:Real};
                       quantity::Symbol=:epsilon, title="Material dispersion", kwargs...)
    x = collect(fs) ./ 1e12
    if quantity === :epsilon
        ε = [𝜀̃(material, f) for f in fs]
        plt = plot(x, real.(ε); label="Re ε", lw=2, xlabel="frequency (THz)",
                   ylabel="permittivity ε", title=title, legend=:best)
        plot!(plt, x, imag.(ε); label="Im ε", lw=2)
    elseif quantity === :index
        nk = [ñ(material, f) for f in fs]
        plt = plot(x, real.(nk); label="n", lw=2, xlabel="frequency (THz)",
                   ylabel="refractive index", title=title, legend=:best)
        plot!(plt, x, imag.(nk); label="k", lw=2)
    else
        error("quantity must be :epsilon or :index, got $quantity")
    end
    plot!(plt; kwargs...)
end


# ---- R / T / A spectra ----
# Plot reflectance, transmittance, and absorptance against an arbitrary x-axis (`xs`, e.g.
# wavelength in nm or frequency in THz) for a vector of TMM or RCWA results. `components`
# selects which curves to draw.
function plot_rta(xs::AbstractVector, results::AbstractVector;
                  xlabel="x", title="Spectrum", components=(:R, :T, :A), kwargs...)
    length(xs) == length(results) ||
        error("xs and results must have equal length ($(length(xs)) vs $(length(results)))")
    R = Float64[]; T = Float64[]; A = Float64[]
    for r in results
        rr, tt, aa = _rta(r)
        push!(R, rr); push!(T, tt); push!(A, aa)
    end
    plt = plot(; xlabel=xlabel, ylabel="efficiency", title=title, legend=:best, kwargs...)
    :R in components && plot!(plt, xs, R; label="R", lw=2)
    :T in components && plot!(plt, xs, T; label="T", lw=2)
    :A in components && plot!(plt, xs, A; label="A", lw=2)
    plt
end


# ---- Reciprocal-space Fourier module ----
# Scatter the retained plane-wave harmonics of a single-lattice RCWA grid in (gx, gy) space,
# highlighting the zeroth order.
function plot_reciprocal_basis(grid::RCWAGrid; title="Reciprocal basis", kwargs...)
    plt = scatter(grid.gx, grid.gy; markersize=5, label="harmonics", aspect_ratio=:equal,
                  xlabel="gx (rad/m)", ylabel="gy (rad/m)", title=title, legend=:best, kwargs...)
    scatter!(plt, [grid.gx[grid.center]], [grid.gy[grid.center]];
             markersize=8, color=:red, label="(0,0)")
    plt
end

# Scatter the joint reciprocal basis G₁ ⊕ … ⊕ G_N of a twisted multilayer stack (the Fourier
# module). With `show_generators=true`, each layer's (rotated) generating vectors b1, b2 are
# drawn as arrows from the origin, making the per-layer sublattices visible.
function plot_reciprocal_basis(grid::TwistedGrid;
                               title="Joint reciprocal basis (Fourier module)",
                               show_generators::Bool=true, kwargs...)
    plt = scatter(grid.gx, grid.gy; markersize=5, label="joint harmonics", aspect_ratio=:equal,
                  xlabel="gx (rad/m)", ylabel="gy (rad/m)", title=title, legend=:best, kwargs...)
    scatter!(plt, [grid.gx[grid.center]], [grid.gy[grid.center]];
             markersize=8, color=:red, label="origin")
    if show_generators
        for i in 1:ndevicelayers(grid)
            c = cos(grid.twist[i]); s = sin(grid.twist[i])
            for (k, b) in enumerate((grid.b1[i], grid.b2[i]))
                vx = c * b[1] - s * b[2]
                vy = s * b[1] + c * b[2]
                plot!(plt, [0.0, vx], [0.0, vy]; lw=2,
                      label=(k == 1 ? "layer $i" : ""))
            end
        end
    end
    plt
end


# ---- Real-space permittivity patterns ----
# Reconstruct and image ε(x, y) = Σ_g ε̂(g, f)·exp(i g·r) of a Fourier periodic pattern. In-plane
# periods are inferred from the smallest non-zero reciprocal components; `n_periods` tiles the
# unit cell, `part` selects which scalar component of ε to show.
function plot_pattern(pattern::FourierPeriodicPattern, f::Real;
                      n_periods::Real=1, npts::Integer=240, part::Symbol=:real,
                      title="ε(x, y)", kwargs...)
    gs = pattern.g_vectors
    coeffs = ε̂(pattern, f)
    gxs = [abs(g[1]) for g in gs if abs(g[1]) > 1e-12]
    gys = [abs(g[2]) for g in gs if abs(g[2]) > 1e-12]
    Lx = isempty(gxs) ? 1.0 : 2π / minimum(gxs)
    Ly = isempty(gys) ? 1.0 : 2π / minimum(gys)
    xs = range(0, n_periods * Lx; length=npts)
    ys = range(0, n_periods * Ly; length=npts)
    Z = Matrix{Float64}(undef, length(ys), length(xs))
    for (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        val = zero(Complex{Float64})
        for k in eachindex(coeffs)
            val += coeffs[k] * cis(gs[k][1] * x + gs[k][2] * y)
        end
        Z[iy, ix] = _partval(val, part)
    end
    heatmap(xs, ys, Z; aspect_ratio=:equal, xlabel="x (m)", ylabel="y (m)",
            title=title, color=:viridis, colorbar_title="ε", kwargs...)
end

plot_pattern(layer::PeriodicLayer, f::Real; kwargs...) = plot_pattern(layer.pattern, f; kwargs...)

# Image a functional (periodic or aperiodic) pattern by sampling ε(x, y, f) over `extent`
# = (xmin, xmax, ymin, ymax).
function plot_pattern(pattern::Union{FunctionalPeriodicPattern,FunctionalAperiodicPattern}, f::Real;
                      extent=(0.0, 1.0, 0.0, 1.0), npts::Integer=240, part::Symbol=:real,
                      title="ε(x, y)", kwargs...)
    xmin, xmax, ymin, ymax = extent
    xs = range(xmin, xmax; length=npts)
    ys = range(ymin, ymax; length=npts)
    Z = [_partval(ε̃_inplane(pattern, x, y, f), part) for y in ys, x in xs]
    heatmap(xs, ys, Z; aspect_ratio=:equal, xlabel="x (m)", ylabel="y (m)",
            title=title, color=:viridis, colorbar_title="ε", kwargs...)
end

# Image a grid-sampled aperiodic pattern directly from its stored data.
function plot_pattern(pattern::GridAperiodicPattern, f::Real=0.0;
                      part::Symbol=:real, title="ε(x, y)", kwargs...)
    Z = [_partval(pattern.ε_data[ix, iy], part)
         for iy in eachindex(pattern.y), ix in eachindex(pattern.x)]
    heatmap(pattern.x, pattern.y, Z; aspect_ratio=:equal, xlabel="x (m)", ylabel="y (m)",
            title=title, color=:viridis, colorbar_title="ε", kwargs...)
end

plot_pattern(layer::AperiodicLayer, f::Real=0.0; kwargs...) = plot_pattern(layer.pattern, f; kwargs...)


# ---- Diffraction orders ----
# Bar charts of the per-order reflectance and transmittance of an RCWA result. Orders whose
# R and T are both below `atol` (typically evanescent) are dropped for readability.
function plot_diffraction_orders(result::RCWAResult; atol::Real=1e-9,
                                 title="Diffraction orders", kwargs...)
    idx = findall(i -> result.R[i] > atol || result.T[i] > atol, eachindex(result.orders))
    isempty(idx) && (idx = collect(eachindex(result.orders)))
    labels = ["($(result.orders[i][1]),$(result.orders[i][2]))" for i in idx]
    xpos = collect(1:length(idx))
    pR = bar(xpos, result.R[idx]; xticks=(xpos, labels), xrotation=60, legend=false,
             ylabel="R", title="reflection", color=:steelblue)
    pT = bar(xpos, result.T[idx]; xticks=(xpos, labels), xrotation=60, legend=false,
             ylabel="T", title="transmission", color=:indianred)
    plot(pR, pT; layout=(1, 2), plot_title=title, kwargs...)
end


# ---- Layer-stack cross-section ----
# Schematic of a layer stack along z: each layer is drawn as a horizontal band whose height is
# proportional to its thickness and whose colour encodes the (real) refractive index at `f`.
# The semi-infinite incidence/transmission half-spaces are drawn with a fixed display height
# (`halfspace_frac` of the finite-stack thickness).
function plot_stack(layers::AbstractVector{<:AbstractLayer}, f::Real;
                    halfspace_frac::Real=0.18, title="Layer stack", kwargs...)
    ns = [real(ñ(material(l), f)) for l in layers]
    finite_t = [thickness(l) for l in layers if isfinite(thickness(l))]
    total = isempty(finite_t) ? 1.0 : sum(finite_t)
    band = halfspace_frac * total
    nmin, nmax = extrema(ns)
    cg = cgrad(:viridis)
    colorof(n) = nmax > nmin ? cg[(n - nmin) / (nmax - nmin)] : cg[0.5]

    plt = plot(; title=title, legend=false, yflip=true, xlims=(0, 1),
               xticks=false, xlabel="", ylabel="z (m)", kwargs...)
    z = 0.0
    for (i, l) in enumerate(layers)
        h = isfinite(thickness(l)) ? thickness(l) : band
        plot!(plt, Shape([0.0, 1.0, 1.0, 0.0], [z, z, z + h, z + h]);
              fillcolor=colorof(ns[i]), linecolor=:black, lw=0.6, label="")
        tag = isfinite(thickness(l)) ? "n=$(round(ns[i], digits=2))" :
              (i == 1 ? "incidence  n=$(round(ns[i], digits=2))" :
                        "substrate  n=$(round(ns[i], digits=2))")
        annotate!(plt, 0.5, z + h / 2, text(tag, 8, :white))
        z += h
    end
    plt
end


# ---- Angle-resolved spectrum / band map ----
# Heatmap of a swept quantity Z over a 2D parameter grid (e.g. angle vs frequency). Z must be
# size (length(ys), length(xs)); rows index ys (frequency axis), columns index xs (the swept
# in-plane parameter). Resonance features trace the photonic bands.
function plot_band_map(xs::AbstractVector, ys::AbstractVector, Z::AbstractMatrix;
                       xlabel="angle θ (deg)", ylabel="frequency (THz)",
                       title="Angle-resolved spectrum", clabel="T", kwargs...)
    size(Z) == (length(ys), length(xs)) ||
        error("Z must be size (length(ys), length(xs)) = ($(length(ys)), $(length(xs))), got $(size(Z))")
    heatmap(xs, ys, Z; xlabel=xlabel, ylabel=ylabel, title=title,
            color=:inferno, colorbar_title=clabel, kwargs...)
end


# ---- Band structures ----
# Plot a (complex) band matrix `B` of size (length(xs), nbands) against the path coordinate
# `xs`. `part` selects the resonance frequency (:real), linewidth (:imag), or magnitude (:abs).
function plot_bands(xs::AbstractVector, B::AbstractMatrix; part::Symbol=:real,
                    xlabel="k", title="Band structure", kwargs...)
    length(xs) == size(B, 1) ||
        error("length(xs) must equal size(B, 1) ($(length(xs)) vs $(size(B, 1)))")
    ylab = part === :real ? "Re(ω)" : part === :imag ? "Im(ω)" : "|ω|"
    plt = plot(; xlabel=xlabel, ylabel=ylab, title=title, legend=false, kwargs...)
    for b in 1:size(B, 2)
        plot!(plt, xs, [_partval(B[p, b], part) for p in 1:size(B, 1)]; lw=2)
    end
    plt
end
