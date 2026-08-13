using LinearAlgebra

# Layer.jl is included (not `using Layer`) to match the flat layout of this repo; it in turn
# loads Material.jl. Load it before TMM references AbstractLayer, material, thickness, and 𝜀̃.
if !isdefined(@__MODULE__, :AbstractLayer)
    include(joinpath(@__DIR__, "Layer.jl"))
end

# Speed of light in vacuum (m/s), matching the constant used in Material.jl. Frequencies f
# are in Hz and lengths (thicknesses) are in meters, so all phases kz·d are dimensionless.
if !isdefined(@__MODULE__, :c₀)
    const c₀ = 299792458.0
end


# ---- Polarization Definition ----
abstract type AbstractPolarization end

struct TE <: AbstractPolarization
    """
    TE()
    Transverse-electric (s) polarization: the electric field is perpendicular to the plane of incidence.
    """
end

struct TM <: AbstractPolarization
    """
    TM()
    Transverse-magnetic (p) polarization: the magnetic field is perpendicular to the plane of incidence.
    """
end


# ---- General accessors for the transfer-matrix method ----
# Bulk relative permittivity seen by the TMM. The method treats every layer as laterally
# homogeneous, so the layer's bulk material model is used. Patterned layers (PeriodicLayer /
# AperiodicLayer without a g = 0 average) intentionally error here: plain TMM is only valid
# for unpatterned stacks.
function 𝜀̃_tmm(layer::AbstractLayer, f::Real)
    𝜀̃(material(layer), f)
end

# Out-of-plane wavevector component kz in a medium of permittivity ε, given the vacuum
# wavenumber k0 and the conserved in-plane wavevector kpar. The branch is chosen so that the
# wave either decays (Im kz ≥ 0) or, when lossless, propagates forward (Re kz ≥ 0) along +z.
function _kz(ε, k0::Real, kpar)
    kz = sqrt(Complex{Float64}(k0^2 * ε - kpar^2))
    if imag(kz) < 0 || (imag(kz) == 0 && real(kz) < 0)
        kz = -kz
    end
    kz
end


# ---- Interface and propagation transfer matrices ----
# 2×2 interface matrix D_{ij} from the Fresnel amplitude coefficients, relating the
# forward/backward field amplitudes on the i-side to those on the j-side:
#     [E⁺_i; E⁻_i] = D_{ij} [E⁺_j; E⁻_j],   D_{ij} = (1/t_{ij}) · [1 r_{ij}; r_{ij} 1].
# εi, εj and kzi, kzj are the permittivities and out-of-plane wavevectors of the two media.
function _interface_matrix(::TE, εi, εj, kzi, kzj)
    r = (kzi - kzj) / (kzi + kzj)
    t = 2 * kzi / (kzi + kzj)
    (1 / t) * [one(r) r; r one(r)]
end

function _interface_matrix(::TM, εi, εj, kzi, kzj)
    r = (εj * kzi - εi * kzj) / (εj * kzi + εi * kzj)
    t = 2 * sqrt(Complex{Float64}(εi * εj)) * kzi / (εj * kzi + εi * kzj)
    (1 / t) * [one(r) r; r one(r)]
end

# 2×2 propagation matrix through a homogeneous layer of thickness d with out-of-plane
# wavevector kz: P = diag(exp(-i·kz·d), exp(+i·kz·d)).
function _propagation_matrix(kz, d::Real)
    δ = kz * d
    Complex{Float64}[exp(-1im * δ) 0; 0 exp(1im * δ)]
end


# ---- TMMResult Definition ----
struct TMMResult
    """
    TMMResult(r, t, R, T, A)
    Result of a transfer-matrix calculation for a single (f, θ, polarization).
    - r: Complex amplitude reflection coefficient
    - t: Complex amplitude transmission coefficient
    - R: Reflectance (power), |r|²
    - T: Transmittance (power), corrected by the impedance ratio of the two half-spaces
    - A: Absorptance, 1 − R − T
    """
    r::Complex{Float64}
    t::Complex{Float64}
    R::Float64
    T::Float64
    A::Float64
end

# Power transmittance from the amplitude coefficient t, including the ratio of longitudinal
# admittances between the transmission and incidence half-spaces.
function _transmittance(::TE, t, εi, εf, kzi, kzf)
    abs2(t) * real(kzf) / real(kzi)
end

function _transmittance(::TM, t, εi, εf, kzi, kzf)
    ni = sqrt(Complex{Float64}(εi))
    nf = sqrt(Complex{Float64}(εf))
    num = real(nf * conj(kzf) / conj(nf))
    den = real(ni * conj(kzi) / conj(ni))
    abs2(t) * num / den
end


# ---- Transfer-matrix solver ----
# The transfer-matrix method for a stack of laterally homogeneous layers under plane-wave
# illumination. The first and last entries of `layers` are the (semi-infinite) incidence and
# transmission half-spaces; their thicknesses are ignored. Inner layers must have finite
# thickness. θ is the angle of incidence (radians) measured from the +z axis in the incidence
# medium, and `pol` selects TE (s) or TM (p) polarization.
function solve(layers::AbstractVector{<:AbstractLayer}, f::Real;
               θ::Real=0.0, pol::AbstractPolarization=TE())
    N = length(layers)
    N >= 2 || error("need at least two layers (incidence and transmission half-spaces), got $N")
    k0 = 2π * f / c₀
    εs = [𝜀̃_tmm(layer, f) for layer in layers]
    kpar = k0 * sqrt(Complex{Float64}(εs[1])) * sin(θ)
    kzs = [_kz(ε, k0, kpar) for ε in εs]

    M = Matrix{Complex{Float64}}(I, 2, 2)
    for i in 1:(N - 1)
        M = M * _interface_matrix(pol, εs[i], εs[i + 1], kzs[i], kzs[i + 1])
        if i + 1 < N
            d = thickness(layers[i + 1])
            isfinite(d) || error("inner layer $(i + 1) must have finite thickness, got $d")
            M = M * _propagation_matrix(kzs[i + 1], d)
        end
    end

    r = M[2, 1] / M[1, 1]
    t = 1 / M[1, 1]
    R = abs2(r)
    T = _transmittance(pol, t, εs[1], εs[N], kzs[1], kzs[N])
    TMMResult(Complex{Float64}(r), Complex{Float64}(t), Float64(R), Float64(T), Float64(1 - R - T))
end


# ---- Convenience wrappers ----
# Scalar reflectance R for the given stack and incidence conditions.
function reflectance(layers::AbstractVector{<:AbstractLayer}, f::Real;
                     θ::Real=0.0, pol::AbstractPolarization=TE())
    solve(layers, f; θ=θ, pol=pol).R
end

# Scalar transmittance T for the given stack and incidence conditions.
function transmittance(layers::AbstractVector{<:AbstractLayer}, f::Real;
                       θ::Real=0.0, pol::AbstractPolarization=TE())
    solve(layers, f; θ=θ, pol=pol).T
end

# Scalar absorptance A = 1 − R − T for the given stack and incidence conditions.
function absorptance(layers::AbstractVector{<:AbstractLayer}, f::Real;
                     θ::Real=0.0, pol::AbstractPolarization=TE())
    solve(layers, f; θ=θ, pol=pol).A
end
