# metasurface_ir — reusable library for the LWIR (8–14 µm) metasurface RCWA workflow.
#
# Include after `using GESTALT`. Provides:
#   • dispersive materials: VO2 (insulating / metallic) and an aluminium mirror;
#   • grayscale-image → thickness map → RCWA layer stack (binary + staircase slicing);
#   • 8–14 µm R/T/A spectral sweep (parallel, TE/TM-averaged by default);
#   • binary shape generators (rectangle / circle / triangle / cross) for the top pattern.
#
# Input convention: a square grayscale image `img` (N×N, values 0–255) encodes the top
# pattern-layer thickness via  gray → (gray/255)·max_t, with pixels ≤ t_min_gray treated as
# air. Binary images (0 = air, one value v > 50 = VO2) are the default input; a grayscale
# (multi-height) image is handled by staircase slicing into uniform-thickness sub-layers.

using DelimitedFiles
using Interpolations

# Directory holding shared material measurement data (measured VO2 n,k tables).
const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))

# ---- Materials -----------------------------------------------------------------------------

# Insulating (monoclinic) VO2: transparent dielectric in the LWIR (n ≈ 2.9 ⇒ ε ≈ 8.4),
# modelled as a single Lorentz oscillator at 17.6 THz (tail of the strong phonons at
# ~17/20 µm).  ε(ω) = ε∞ + Δε·ω₀² / (ω₀² − ω² − i·γ₀·ω).
vo2_insulating_material() = LorentzMaterial(8.4, [(2.0, 2π * 17.6e12, 2π * 1.6e12)])

# Metallic (rutile) VO2: a lossy "bad metal" in the LWIR, Drude model
#   ε(ω) = ε∞ − ωₚ² / (ω² + i·γₚ·ω).
vo2_metallic_material() = DrudeMaterial(1.51e15, 5.75e14, 9.0)

# Aluminium mirror: Drude metal (ωₚ ≈ 2.4e16 rad/s, γ ≈ 1.24e14 rad/s). At the default 200 nm
# thickness it is ~10 skin depths thick in the LWIR, i.e. effectively opaque (T ≈ 0).
aluminum_material() = DrudeMaterial(2.4e16, 1.24e14, 1.0)


# ---- Sampled (measured) VO2 from CSV -------------------------------------------------------
# Measured VO2 refractive-index tables (n,k vs wavelength). CSV columns are
# `wavelengthmicrons, Reindex, Imindex` in descending wavelength. `SampledDataMaterial`
# stores wavelength in nm ascending, so the loader converts units and reverses the order.
# VO2-I = insulating (monoclinic), VO2-M = metallic (rutile), over 0.4–30 µm.

# Generic loader for a three-column (λ, n, k) CSV with a header row.
#   • λ_unit: unit of the wavelength column in metres (default 1e-6 ⇒ µm).
#   • interp_kind: interpolation scheme (default linear — robust for measured n,k).
function load_sampled_nk_csv(path::String; λ_unit::Real=1e-6,
                             interp_kind=Gridded(Linear()))
    raw = readdlm(path, ',', skipstart=1)
    size(raw, 2) == 3 ||
        error("CSV must have 3 columns (wavelength, n, k); got $(size(raw, 2))")
    λ_nm = Float64.(raw[:, 1]) .* (λ_unit * 1e9)      # → nm
    n = Float64.(raw[:, 2])
    k = Float64.(raw[:, 3])
    if λ_nm[1] > λ_nm[end]                            # descending ⇒ reverse to ascending
        λ_nm = reverse(λ_nm); n = reverse(n); k = reverse(k)
    end
    keep = [i == 1 || λ_nm[i] > λ_nm[i - 1] for i in eachindex(λ_nm)]
    if !all(keep)                                     # drop duplicate wavelength points
        λ_nm = λ_nm[keep]; n = n[keep]; k = k[keep]
    end
    SampledDataMaterial(λ_nm, n, k; interp_kind=interp_kind)
end

# Insulating VO2 from the measured table (data/VO2-I.csv).
vo2_insulating_sampled() = load_sampled_nk_csv(joinpath(_DATA_DIR, "VO2-I.csv"))

# Metallic VO2 from the measured table (data/VO2-M.csv).
vo2_metallic_sampled() = load_sampled_nk_csv(joinpath(_DATA_DIR, "VO2-M.csv"))


# ---- Grayscale → thickness ----------------------------------------------------------------

# Map one grayscale value to a material thickness: gray ∈ [0, max_gray] → [0, max_t]. Pixels at
# or below `t_min_gray` are air (thickness 0). `max_gray` is the grayscale level that maps to
# `max_t` (default 255, the 8-bit convention).
gray_to_thickness(gray::Real, max_t::Real=4e-6, t_min_gray::Real=50, max_gray::Real=255) =
    gray > t_min_gray ? (gray / max_gray) * max_t : 0.0

# Full thickness map (meters) of a square grayscale image.
function thickness_profile(img::AbstractMatrix, max_t::Real=4e-6, t_min_gray::Real=50,
                           max_gray::Real=255)
    [gray_to_thickness(img[i, j], max_t, t_min_gray, max_gray)
     for i in 1:size(img, 1), j in 1:size(img, 2)]
end


# ---- Binary image → Fourier pattern -------------------------------------------------------

# Fourier amplitudes of a binary in-plane pattern given directly as a square indicator matrix
# (1 = material, 0 = background) on the unit cell of a square lattice with period `p`:
#   ε(r, f) = ε_out(f) + (ε_in(f) − ε_out(f))·indicator(r),
#   ε̂(m,n)  = ε_out·δ_{m0,n0} + (ε_in − ε_out)·F(m,n),
#   F(m,n)  = (1/N²) Σ_{i,j} ind[i,j]·exp(−i2π(m(i−1)+n(j−1))/N).
# The stored g-vectors are exact multiples of the reciprocal basis (m·2π/p, n·2π/p), so the
# RCWA convolution matrix can look them up by integer order differences. `mmax`, `nmax` should
# be ≥ 2× the RCWA truncation so every convolution difference g_p − g_q is available.
function fourier_pattern_from_indicator(ind::AbstractMatrix, p::Real,
                                        mat_in::AbstractMaterial, mat_out::AbstractMaterial;
                                        mmax::Integer, nmax::Integer)
    N = size(ind, 1)
    size(ind, 2) == N || error("indicator must be square (N×N), got $(size(ind))")
    b = 2π / p
    gvs = NTuple{2,Float64}[]
    coeffs = Any[]
    for n in (-nmax):nmax, m in (-mmax):mmax
        acc = 0.0 + 0.0im
        @inbounds for j in 1:N, i in 1:N
            acc += Float64(ind[i, j]) * cis(-2π * (m * (i - 1) / N + n * (j - 1) / N))
        end
        F = acc / N^2
        is0 = (m == 0 && n == 0)
        push!(gvs, (m * b, n * b))
        push!(coeffs, f -> (is0 ? 𝜀̃(mat_out, f) : 0.0 + 0.0im) +
                           (𝜀̃(mat_in, f) - 𝜀̃(mat_out, f)) * F)
    end
    FourierPeriodicPattern(gvs, coeffs)
end


# ---- Image → RCWA layer stack -------------------------------------------------------------

# Build the patterned top layer(s) from a grayscale image. A binary image (single non-zero
# height) gives one `PeriodicLayer` of that thickness; a multi-height image is approximated by
# `nslices` staircase sub-layers of equal thickness, each a binary VO2/air pattern.
function build_pattern_layers(img::AbstractMatrix, mat_vo2::AbstractMaterial,
                              mat_bg::AbstractMaterial=NKMaterial(1.0);
                              p::Real=5e-6, max_t::Real=4e-6, t_min_gray::Real=50,
                              max_gray::Real=255, N::Integer=8, nslices::Integer=8)
    h = thickness_profile(img, max_t, t_min_gray, max_gray)
    heights = sort!(unique(filter(>(0.0), h)))
    isempty(heights) && return AbstractLayer[]          # all air ⇒ no top layer
    mm = 2N
    if length(heights) == 1                              # binary image: single layer
        ind = Float64.(h .> 0.0)
        pat = fourier_pattern_from_indicator(ind, p, mat_vo2, mat_bg; mmax=mm, nmax=mm)
        return AbstractLayer[PeriodicLayer(mat_vo2, heights[1], pat)]
    end
    # grayscale: staircase slicing into equal-thickness binary sub-layers
    hmax = heights[end]
    Δz = hmax / nslices
    layers = AbstractLayer[]
    for k in 1:nslices
        zbot = (k - 1) * Δz
        ind = Float64.(h .>= zbot)
        any(>(0.0), ind) || continue
        pat = fourier_pattern_from_indicator(ind, p, mat_vo2, mat_bg; mmax=mm, nmax=mm)
        push!(layers, PeriodicLayer(mat_vo2, Δz, pat))
    end
    layers
end

# Full metasurface stack: air | VO2 pattern | (optional uniform dielectric spacers...) |
# aluminium mirror (200 nm) | air. The `spacers` interface is reserved for future multi-layer
# designs: a list of `HomogeneousLayer`s inserted between the pattern and the mirror.
# Returns `(grid, layers)` where `grid` is the RCWA truncation grid for the pattern.
function build_stack(img::AbstractMatrix, mat_vo2::AbstractMaterial,
                     mat_al::AbstractMaterial=aluminum_material();
                     p::Real=5e-6, max_t::Real=4e-6, t_min_gray::Real=50,
                     N::Integer=8, nslices::Integer=8,
                     spacers::AbstractVector{<:AbstractLayer}=AbstractLayer[],
                     d_al::Real=200e-9)
    top = build_pattern_layers(img, mat_vo2; p=p, max_t=max_t, t_min_gray=t_min_gray,
                               N=N, nslices=nslices)
    air = HomogeneousLayer(NKMaterial(1.0), Inf)
    al  = HomogeneousLayer(mat_al, d_al)
    layers = AbstractLayer[air, top..., spacers..., al, air]
    grid = RCWAGrid((2π / p, 0.0), (0.0, 2π / p); Nx=N, Ny=N)
    (grid, layers)
end


# ---- Spectral sweep -----------------------------------------------------------------------

# R/T/A spectra over the wavelengths `λs` (meters) at fixed incidence. Parallel over λ;
# call `BLAS.set_num_threads(1)` in the driver to avoid oversubscription. With
# `unpolarized=true` the TE and TM results are averaged (natural/unpolarised light).
function sweep_spectrum(layers::AbstractVector{<:AbstractLayer}, grid::RCWAGrid,
                        λs::AbstractVector{<:Real}; θ::Real=0.0, φ::Real=0.0,
                        unpolarized::Bool=true)
    nλ = length(λs)
    R = zeros(nλ); T = zeros(nλ); A = zeros(nλ)
    Threads.@threads for k in eachindex(λs)
        f = c₀ / λs[k]
        if unpolarized
            rTE = solve(layers, grid, f; θ=θ, φ=φ, pol_TE=1.0, pol_TM=0.0)
            rTM = solve(layers, grid, f; θ=θ, φ=φ, pol_TE=0.0, pol_TM=1.0)
            R[k] = 0.5 * (rTE.Rtot + rTM.Rtot)
            T[k] = 0.5 * (rTE.Ttot + rTM.Ttot)
            A[k] = 0.5 * (rTE.A + rTM.A)
        else
            r = solve(layers, grid, f; θ=θ, φ=φ, pol_TE=1.0, pol_TM=0.0)
            R[k] = r.Rtot; T[k] = r.Ttot; A[k] = r.A
        end
    end
    (R, T, A)
end

# Detailed polarized sweep: returns separate TE and TM spectra so the caller can either
# average them (unpolarized light) or return a single polarization. Parallel over λ.
function sweep_spectrum_detailed(layers::AbstractVector{<:AbstractLayer}, grid::RCWAGrid,
                                 λs::AbstractVector{<:Real}; θ::Real=0.0, φ::Real=0.0)
    nλ = length(λs)
    RTE = zeros(nλ); TTE = zeros(nλ); ATE = zeros(nλ)
    RTM = zeros(nλ); TTM = zeros(nλ); ATM = zeros(nλ)
    Threads.@threads for k in eachindex(λs)
        f = c₀ / λs[k]
        rTE = solve(layers, grid, f; θ=θ, φ=φ, pol_TE=1.0, pol_TM=0.0)
        rTM = solve(layers, grid, f; θ=θ, φ=φ, pol_TE=0.0, pol_TM=1.0)
        RTE[k] = rTE.Rtot; TTE[k] = rTE.Ttot; ATE[k] = rTE.A
        RTM[k] = rTM.Rtot; TTM[k] = rTM.Ttot; ATM[k] = rTM.A
    end
    (R_TE=RTE, T_TE=TTE, A_TE=ATE, R_TM=RTM, T_TM=TTM, A_TM=ATM)
end


# ---- Binary shape generators (N×N grayscale, values 0 / `gray`) ---------------------------

# Point-in-triangle test in fractional cell coordinates (all 2-tuples in [0,1)²) via the
# sign-of-cross-product (barycentric) test; the boundary counts as inside.
function _in_triangle(p, a, b, c)
    d1 = (p[1] - b[1]) * (a[2] - b[2]) - (a[1] - b[1]) * (p[2] - b[2])
    d2 = (p[1] - c[1]) * (b[2] - c[2]) - (b[1] - c[1]) * (p[2] - c[2])
    d3 = (p[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (p[2] - a[2])
    !((d1 < 0 || d2 < 0 || d3 < 0) && (d1 > 0 || d2 > 0 || d3 > 0))
end

# Blank (all-air) image.
blank_image(N::Integer=64) = zeros(UInt8, N, N)

# Axis-aligned rectangle; (x0, y0) lower-left corner, size (w, h), in fractional coords.
function rect_image(x0::Real, y0::Real, w::Real, h::Real; gray::Integer=0xff, N::Integer=64)
    img = zeros(UInt8, N, N)
    for j in 1:N, i in 1:N
        x = (i - 0.5) / N; y = (j - 0.5) / N
        (x0 <= x < x0 + w && y0 <= y < y0 + h) && (img[i, j] = gray)
    end
    img
end

# Filled circle; (x0, y0) centre, r radius, in fractional coords (not wrapped across the cell).
function circle_image(x0::Real, y0::Real, r::Real; gray::Integer=0xff, N::Integer=64)
    img = zeros(UInt8, N, N)
    for j in 1:N, i in 1:N
        x = (i - 0.5) / N; y = (j - 0.5) / N
        ((x - x0)^2 + (y - y0)^2 <= r^2) && (img[i, j] = gray)
    end
    img
end

# Filled triangle from three vertices (fractional coords).
function triangle_image(a, b, c; gray::Integer=0xff, N::Integer=64)
    img = zeros(UInt8, N, N)
    for j in 1:N, i in 1:N
        p = ((i - 0.5) / N, (j - 0.5) / N)
        _in_triangle(p, a, b, c) && (img[i, j] = gray)
    end
    img
end

# Plus/cross shape: a vertical and a horizontal bar of width `arm` inside the (x0,y0,w,h) box.
function cross_image(x0::Real, y0::Real, w::Real, h::Real, arm::Real;
                     gray::Integer=0xff, N::Integer=64)
    img = zeros(UInt8, N, N)
    for j in 1:N, i in 1:N
        x = (i - 0.5) / N; y = (j - 0.5) / N
        vertical   = abs(x - (x0 + w / 2)) <= arm / 2 && y0 <= y < y0 + h
        horizontal = abs(y - (y0 + h / 2)) <= arm / 2 && x0 <= x < x0 + w
        (vertical || horizontal) && (img[i, j] = gray)
    end
    img
end

# Named-shape convenience: returns a binary grayscale image (0 / gray) for one of
# "rectangle", "circle", "triangle" or "cross".
function shape_image(shape::String; N::Integer=64, gray::Integer=0xff)
    s = lowercase(strip(shape))
    if s in ("rectangle", "rect")
        return rect_image(0.15, 0.25, 0.60, 0.40; gray=gray, N=N)
    elseif s in ("circle", "disk")
        return circle_image(0.50, 0.50, 0.34; gray=gray, N=N)
    elseif s in ("triangle", "tri")
        return triangle_image((0.12, 0.12), (0.88, 0.18), (0.42, 0.82); gray=gray, N=N)
    elseif s in ("cross", "plus")
        return cross_image(0.20, 0.20, 0.60, 0.60, 0.22; gray=gray, N=N)
    else
        error("unknown shape '$shape'; use rectangle / circle / triangle / cross")
    end
end
