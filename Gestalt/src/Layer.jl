using Interpolations


if !isdefined(@__MODULE__, :AbstractMaterial)
    include(joinpath(@__DIR__, "Material.jl"))
end


# ---- AbstractLayer Definition ----
abstract type AbstractLayer end

# ---- Base layer properties (material + thickness) ----
function material(layer::AbstractLayer)
    error("material(::$(typeof(layer))) not implemented")
end

function thickness(layer::AbstractLayer)
    error("thickness(::$(typeof(layer))) not implemented")
end

function 𝜀̃_material(layer::AbstractLayer, f::Real)
    𝜀̃(material(layer), f)
end

function ñ_material(layer::AbstractLayer, f::Real)
    ñ(material(layer), f)
end


# ---- Homogeneous layers ----
abstract type AbstractHomogeneousLayer <: AbstractLayer end

struct HomogeneousLayer <: AbstractHomogeneousLayer
    """
    HomogeneousLayer(material::AbstractMaterial, thickness::Real)
    Laterally uniform layer: ε(r) = ε̃(material, f) in the layer plane.
    - material: Bulk material model (see `AbstractMaterial` in `Material.jl`)
    - thickness: Layer thickness along z (must be non-negative; may be `Inf` for a semi-infinite half-space)
    """
    material::AbstractMaterial
    thickness::Float64
    function HomogeneousLayer(material::AbstractMaterial, thickness::Real)
        thickness >= 0 || error("thickness must be non-negative, got $thickness")
        new(material, Float64(thickness))
    end
end

material(layer::HomogeneousLayer) = layer.material
thickness(layer::HomogeneousLayer) = layer.thickness

function 𝜀̃(layer::HomogeneousLayer, f::Real)
    𝜀̃(material(layer), f)
end

function ñ(layer::HomogeneousLayer, f::Real)
    ñ(material(layer), f)
end


# ---- Periodic in-plane patterns ----
abstract type AbstractPeriodicPattern end

function ε̂(pattern::AbstractPeriodicPattern, f::Real)
    error("ε̂(::$(typeof(pattern)), f) not implemented")
end

function g_vectors(pattern::AbstractPeriodicPattern)
    error("g_vectors(::$(typeof(pattern))) not implemented")
end

struct FourierPeriodicPattern <: AbstractPeriodicPattern
    """
    FourierPeriodicPattern(g_vectors, ε̂_coefficients)
    RCWA-style periodic dielectric: ε(r) = Σ_g ε̂(g, f) · exp(i g·r).
    - g_vectors: Reciprocal lattice vectors (gx, gy) in the layer plane
    - ε̂_coefficients: Fourier amplitudes ε̂(g, f); each entry is a `Complex` constant or a `Function` f -> Complex
    Note:
        * The zeroth harmonic (g = 0) gives the spatial average when present.
        * Constituent dispersion is often encoded in the ε̂(g, f) functions; `material` on the parent `PeriodicLayer` is the reference/host material model.
    """
    g_vectors::Vector{NTuple{2,Float64}}
    ε̂_coefficients::Vector{Any}
    function FourierPeriodicPattern(g_vectors::AbstractVector{<:NTuple{2,<:Real}},
                                    ε̂_coefficients::AbstractVector)
        length(g_vectors) == length(ε̂_coefficients) ||
            error("g_vectors and ε̂_coefficients must have the same length")
        isempty(g_vectors) && error("FourierPeriodicPattern requires at least one (g, ε̂) pair")
        g = [(Float64(gx), Float64(gy)) for (gx, gy) in g_vectors]
        coeffs = Any[]
        for (i, c) in enumerate(ε̂_coefficients)
            if c isa Function
                push!(coeffs, c)
            elseif c isa Complex
                push!(coeffs, Complex{Float64}(c))
            elseif c isa Real
                push!(coeffs, Complex{Float64}(c))
            else
                error("ε̂ coefficient $i must be Complex, Real, or Function, got $(typeof(c))")
            end
        end
        new(g, coeffs)
    end
end

g_vectors(pattern::FourierPeriodicPattern) = pattern.g_vectors

function ε̂(pattern::FourierPeriodicPattern, f::Real)
    [_eval_ε̂_coeff(c, f) for c in pattern.ε̂_coefficients]
end

function _eval_ε̂_coeff(c, f::Real)
    c isa Function ? Complex{Float64}(c(f)) : c::Complex{Float64}
end

struct PixelPeriodicPattern <: AbstractPeriodicPattern
    """
    PixelPeriodicPattern(indicator::AbstractMatrix, period::Real,
                         mat_in::AbstractMaterial, mat_out::AbstractMaterial)
    在正方形像素网格上给出的周期面内图案（1=材料 in，0=材料 out）。用于实空间 FDFD
    面内本征求解路径（`src/FDFD.jl`）：ε 在像素上精确对角、导数用稀疏有限差分，避免
    Fourier RCWA 对细结构的 Gibbs/截断不收敛问题。indicator 为 Npix×Npix 矩阵。
    """
    indicator::Matrix{Float64}
    period::Float64
    mat_in::AbstractMaterial
    mat_out::AbstractMaterial
    function PixelPeriodicPattern(indicator::AbstractMatrix{<:Real}, period::Real,
                                  mat_in::AbstractMaterial, mat_out::AbstractMaterial)
        size(indicator, 1) == size(indicator, 2) ||
            error("indicator must be square (Npix×Npix), got $(size(indicator))")
        period > 0 || error("period must be positive")
        new(Float64.(indicator), Float64(period), mat_in, mat_out)
    end
end

function ε̃_inplane(pattern::PixelPeriodicPattern, x::Real, y::Real, f::Real)
    N = size(pattern.indicator, 1)
    i = clamp(floor(Int, x / pattern.period * N) + 1, 1, N)
    j = clamp(floor(Int, y / pattern.period * N) + 1, 1, N)
    ε_in = 𝜀̃(pattern.mat_in, f)
    ε_out = 𝜀̃(pattern.mat_out, f)
    ε_out + (ε_in - ε_out) * pattern.indicator[i, j]
end

struct FunctionalPeriodicPattern{F<:Function} <: AbstractPeriodicPattern
    """
    FunctionalPeriodicPattern(ε_func::Function)
    Periodic in-plane permittivity profile supplied explicitly.
    - ε_func: (x, y, f) -> Complex permittivity; must be periodic in (x, y)
    Note: Reciprocal vectors are not stored; use `FourierPeriodicPattern` when a Fourier list is required for RCWA.
    """
    ε_func::F
end

function ε̃_inplane(pattern::FunctionalPeriodicPattern, x::Real, y::Real, f::Real)
    Complex{Float64}(pattern.ε_func(x, y, f))
end


# ---- Periodic layers ----
abstract type AbstractPeriodicLayer <: AbstractLayer end

struct PeriodicLayer <: AbstractPeriodicLayer
    """
    PeriodicLayer(material::AbstractMaterial, thickness::Real, pattern::AbstractPeriodicPattern)
    Layer uniform along z with a periodic in-plane dielectric pattern (RCWA-compatible).
    - material: Reference/host material model (see `AbstractMaterial` in `Material.jl`)
    - thickness: Layer thickness along z (must be non-negative)
    - pattern: In-plane periodic permittivity pattern (e.g. `FourierPeriodicPattern`)
    """
    material::AbstractMaterial
    thickness::Float64
    pattern::AbstractPeriodicPattern
    function PeriodicLayer(material::AbstractMaterial,
                           thickness::Real,
                           pattern::AbstractPeriodicPattern)
        thickness >= 0 || error("thickness must be non-negative, got $thickness")
        new(material, Float64(thickness), pattern)
    end
end

material(layer::PeriodicLayer) = layer.material
thickness(layer::PeriodicLayer) = layer.thickness
pattern(layer::PeriodicLayer) = layer.pattern

function ε̂(layer::PeriodicLayer, f::Real)
    ε̂(layer.pattern, f)
end

function g_vectors(layer::PeriodicLayer)
    g_vectors(layer.pattern)
end

function 𝜀̃(layer::PeriodicLayer, f::Real)
    pattern = layer.pattern
    if pattern isa FourierPeriodicPattern
        for (g, c) in zip(pattern.g_vectors, pattern.ε̂_coefficients)
            g == (0.0, 0.0) && return _eval_ε̂_coeff(c, f)
        end
        error("PeriodicLayer has no g = 0 Fourier component; use ε̂(layer, f) or 𝜀̃_material(layer, f)")
    else
        error("𝜀̃(::PeriodicLayer, f) requires a FourierPeriodicPattern with g = 0, or use ε̃_inplane(layer, x, y, f)")
    end
end

function ε̃_inplane(layer::PeriodicLayer, x::Real, y::Real, f::Real)
    layer.pattern isa FunctionalPeriodicPattern ||
        error("ε̃_inplane(::PeriodicLayer, ...) requires FunctionalPeriodicPattern")
    ε̃_inplane(layer.pattern, x, y, f)
end


# ---- Aperiodic in-plane patterns ----
abstract type AbstractAperiodicPattern end

function ε̃_inplane(pattern::AbstractAperiodicPattern, x::Real, y::Real, f::Real)
    error("ε̃_inplane(::$(typeof(pattern)), x, y, f) not implemented")
end

struct FunctionalAperiodicPattern{F<:Function} <: AbstractAperiodicPattern
    """
    FunctionalAperiodicPattern(ε_func::Function)
    Aperiodic in-plane permittivity profile.
    - ε_func: (x, y, f) -> Complex permittivity
    """
    ε_func::F
end

function ε̃_inplane(pattern::FunctionalAperiodicPattern, x::Real, y::Real, f::Real)
    Complex{Float64}(pattern.ε_func(x, y, f))
end

struct GridAperiodicPattern <: AbstractAperiodicPattern
    """
    GridAperiodicPattern(x::AbstractVector{<:Real},
                         y::AbstractVector{<:Real},
                         ε::AbstractMatrix{<:Complex};
                         interp_kind=BSpline(Cubic(Line(OnGrid()))))
    Aperiodic permittivity given on a regular (x, y) grid with bilinear/cubic interpolation.
    - x, y: Grid coordinates (must be strictly increasing)
    - ε: Complex permittivity samples, size length(x) × length(y)
    - interp_kind: Interpolation method (default cubic spline; requires `Interpolations.jl`)
    """
    x::Vector{Float64}
    y::Vector{Float64}
    ε_data::Matrix{Complex{Float64}}
    interp_ε

    function GridAperiodicPattern(x::AbstractVector{<:Real},
                                  y::AbstractVector{<:Real},
                                  ε::AbstractMatrix{<:Complex};
                                  interp_kind=BSpline(Cubic(Line(OnGrid()))))
        length(x) == size(ε, 1) || error("length(x) must match size(ε, 1)")
        length(y) == size(ε, 2) || error("length(y) must match size(ε, 2)")
        x_vals = Float64.(x)
        y_vals = Float64.(y)
        _is_strictly_increasing(x_vals) || error("x grid must be strictly increasing")
        _is_strictly_increasing(y_vals) || error("y grid must be strictly increasing")
        ε_vals = Complex{Float64}.(ε)
        itp = interpolate((x_vals, y_vals), ε_vals, interp_kind)
        new(x_vals, y_vals, ε_vals, scale(itp, x_vals, y_vals))
    end
end

function _is_strictly_increasing(v::AbstractVector{<:Real})
    length(v) < 2 && return true
    all(i -> v[i] < v[i + 1], 1:(length(v) - 1))
end

function ε̃_inplane(pattern::GridAperiodicPattern, x::Real, y::Real, f::Real)
    # f is unused for static grids; kept for a uniform (x, y, f) API
    pattern.interp_ε(x, y)
end


# ---- Aperiodic layers ----
abstract type AbstractAperiodicLayer <: AbstractLayer end

struct AperiodicLayer <: AbstractAperiodicLayer
    """
    AperiodicLayer(material::AbstractMaterial, thickness::Real, pattern::AbstractAperiodicPattern)
    Layer uniform along z with a non-periodic in-plane dielectric pattern (real-space CMT / custom grids).
    - material: Reference material model (see `AbstractMaterial` in `Material.jl`)
    - thickness: Layer thickness along z (must be non-negative)
    - pattern: In-plane aperiodic permittivity pattern (e.g. `FunctionalAperiodicPattern`, `GridAperiodicPattern`)
    """
    material::AbstractMaterial
    thickness::Float64
    pattern::AbstractAperiodicPattern
    function AperiodicLayer(material::AbstractMaterial,
                            thickness::Real,
                            pattern::AbstractAperiodicPattern) 
        thickness >= 0 || error("thickness must be non-negative, got $thickness")
        new(material, Float64(thickness), pattern)
    end
end

material(layer::AperiodicLayer) = layer.material
thickness(layer::AperiodicLayer) = layer.thickness
pattern(layer::AperiodicLayer) = layer.pattern

function ε̃_inplane(layer::AperiodicLayer, x::Real, y::Real, f::Real)
    ε̃_inplane(layer.pattern, x, y, f)
end
