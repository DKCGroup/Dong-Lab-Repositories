using Interpolations
using DelimitedFiles


# ---- AbstractMaterial Definition ----
abstract type AbstractMaterial end

# ---- General Definition of 𝜀̃ and ñ for AbstractMaterial ----
function 𝜀̃(material::AbstractMaterial, f::Real)
    error("𝜀̃(::$(typeof(material)), f) not implemented")
end

function ñ(material::AbstractMaterial, f::Real)
    ε = 𝜀̃(material, f)
    n = sqrt((abs(ε) + real(ε)) / 2)
    k = sqrt((abs(ε) - real(ε)) / 2)
    Complex{Float64}(n, k)
end


# ---- NKMaterial Definition ----
struct NKMaterial <: AbstractMaterial
    """
    NKMaterial(n::Real, k::Real=0.0) 
    Material defined by constant refractive index n and extinction coefficient k.
    - n: Refractive index (must be non-negative)
    - k: Extinction coefficient (must be non-negative)
    """
    n::Float64
    k::Float64
    function NKMaterial(n::Real, k::Real=0.0)
        n >= 0 || error("n must be non-negative, got $n")
        k >= 0 || error("k must be non-negative, got $k")
        new(Float64(n), Float64(k))
    end
end

# For NKMaterial, 𝜀̃ and ñ can be directly computed from n and k without frequency dependence.
𝜀̃(material::NKMaterial, f::Real) = Complex{Float64}(material.n + 1im*material.k)^2
ñ(material::NKMaterial, f::Real) = Complex{Float64}(material.n, material.k)

# ---- FunctionalNKMaterial Definition ----
struct FunctionalNKMaterial{Func<:Function} <: AbstractMaterial
    """
    FunctionalNKMaterial(n_func::Function, k_func::Function)
    Material defined by frequency-dependent refractive index and extinction coefficient.
    - n_func: Function that takes frequency f and returns refractive index n(f) (must be non-negative for all f)
    - k_func: Function that takes frequency f and returns extinction coefficient k(f) (must be non-negative for all f)
    """
    n_func::Func
    k_func::Func
end

# For FunctionalNKMaterial, 𝜀̃ and ñ are computed by evaluating n_func and k_func at the given frequency f.
𝜀̃(material::FunctionalNKMaterial, f::Real) = Complex{Float64}(material.n_func(f) + 1im*material.k_func(f))^2
ñ(material::FunctionalNKMaterial, f::Real) = Complex{Float64}(material.n_func(f), material.k_func(f))

# ---- FunctionalEpsilonMaterial Definition ----
struct FunctionalEpsilonMaterial{Func<:Function} <: AbstractMaterial
    """
    FunctionalEpsilonMaterial(ε_func::Function)
    Material defined directly by a frequency-dependent complex relative permittivity.
    - ε_func: Function f -> ε̃(f) (frequency f in Hz; returns a complex permittivity)
    Use this for dispersive models naturally expressed in ε rather than in (n, k) — e.g.
    multi-oscillator phonon fits with complex oscillator strengths or a complex ε∞ (which the
    real-valued `LorentzMaterial` cannot represent). The refractive index follows from the
    generic `ñ(::AbstractMaterial, f)` fallback (e^{-iωt} convention, Im ε ≥ 0 for loss).
    """
    ε_func::Func
end

𝜀̃(material::FunctionalEpsilonMaterial, f::Real) = Complex{Float64}(material.ε_func(f))


# ---- LorentzMaterial Definition ----
struct LorentzMaterial <: AbstractMaterial
    """
    LorentzMaterial(𝜀̃∞::Real, oscillators::AbstractVector{<:Tuple{Real,Real,Real}})
    Material defined by a Lorentz oscillator model.
    - 𝜀̃∞: High-frequency permittivity (must be positive)
    - oscillators: Vector of tuples (Δ𝜀̃, ω₀, γ₀) defining each Lorentz oscillator
        - Δ𝜀̃: Oscillator strength (must be > 0)
        - ω₀: Resonance frequency (must be > 0)
        - γ₀: Damping coefficient (must be ≥ 0)
    The permittivity is given by:
        𝜀̃(ω) = 𝜀̃∞ + Σ [ Δ𝜀̃ * ω₀^2 / (ω₀^2 - ω^2 - i*γ₀*ω) ]
    where the sum is taken over all oscillators.
    Note: 
        * For a lossless Lorentz material, set γ₀ = 0 for all oscillators.
        * For a perfect electric conductor (PEC), set Δ𝜀̃ → ∞ for at least one oscillator.
    """
    𝜀̃∞::Float64
    oscillators::Vector{Tuple{Float64,Float64,Float64}}
    function LorentzMaterial(𝜀̃∞::Real, oscillators::AbstractVector{<:Tuple{Real,Real,Real}})
        𝜀̃∞ > 0 || error("High-frequency permittivity 𝜀̃∞ must be positive, got $𝜀̃∞")
        for (i, (Δ𝜀̃, ω₀, γ₀)) in enumerate(oscillators)
            Δ𝜀̃ > 0 || error("oscillator $i: Lorentz permittivity Δ𝜀̃ must be > 0, got $Δ𝜀̃")
            ω₀ > 0 || error("oscillator $i: Lorentz Oscillating frequency ω₀ must be > 0, got $ω₀")
            γ₀ >= 0 || error("oscillator $i: LorentzDamping coefficient γ₀ must be ≥ 0, got $γ₀")
        end
        new(Float64(𝜀̃∞), [(Float64(Δ𝜀̃), Float64(ω₀), Float64(γ₀)) for (Δ𝜀̃, ω₀, γ₀) in oscillators])
    end
end

# The permittivity of a LorentzMaterial is computed by summing the contributions from all oscillators according to the formula given in the docstring.
function 𝜀̃(material::LorentzMaterial, f::Real)
    ω = 2π * f
    eps = material.𝜀̃∞ + 0.0im
    for (Δ𝜀̃, ω₀, γ₀) in material.oscillators
        eps += Complex{Float64}(Δ𝜀̃ * ω₀^2 / (ω₀^2 - ω^2 - 1im*γ₀*ω))
    end
    eps
end


# ---- DrudeMaterial Definition ----
struct DrudeMaterial <: AbstractMaterial
    """
    DrudeMaterial(ωₚ::Real, γₚ::Real, 𝜀̃∞::Real=1.0)
    Material defined by the Drude model for free-electron metals.
        - ωₚ: Plasma frequency (must be ≥ 0)
        - γₚ: Damping coefficient (must be ≥ 0)
        - 𝜀̃∞: High-frequency permittivity (must be positive, default is 1.0)
    The permittivity is given by:
        𝜀̃(ω) = 𝜀̃∞ - ωₚ^2 / (ω^2 + i*γₚ*ω)
    where ω = 2πf is the angular frequency.
    Note: 
        * For a lossless Drude metal, set γₚ = 0. 
        * For a perfect electric conductor (PEC), set ωₚ → ∞.
    """
    ωₚ::Float64
    γₚ::Float64
    𝜀̃∞::Float64
    function DrudeMaterial(ωₚ::Real, γₚ::Real, 𝜀̃∞::Real=1.0)
        𝜀̃∞ > 0 || error("Drude permittivity 𝜀̃∞ must be positive, got $𝜀̃∞")
        ωₚ >= 0 || error("Drude pole frequency ωₚ must be ≥0, got $ωₚ")
        γₚ >= 0  || error("Drude pole damping coefficient γₚ must be ≥0, got $γₚ")
        new(Float64(ωₚ), Float64(γₚ), Float64(𝜀̃∞))
    end
end

# The permittivity of a DrudeMaterial is computed using the formula given in the docstring, which accounts for the plasma frequency, damping, and high-frequency permittivity.
function 𝜀̃(material::DrudeMaterial, f::Real)
    ω = 2π * f
    material.𝜀̃∞ - material.ωₚ^2 / (ω^2 + 1im*material.γₚ*ω)
end


# ---- DebyeMaterial Definition ----
struct DebyeMaterial <: AbstractMaterial
    """
    DebyeMaterial(𝜀̃ₛ::Real, 𝜀̃∞::Real, τ::Real)
    Material defined by the Debye relaxation model for polar dielectrics.
        - 𝜀̃ₛ: Static permittivity (must be positive)
        - 𝜀̃∞: High-frequency permittivity (must be positive)
        - τ: Relaxation time (must be positive)
    The permittivity is given by:
        𝜀̃(ω) = 𝜀̃∞ + (𝜀̃ₛ - 𝜀̃∞) / (1 + i*ω*τ)
    where ω = 2πf is the angular frequency.
    Note: 
        * For a lossless Debye material, set τ → 0. 
        * For a perfect electric insulator, set 𝜀̃ₛ → ∞.
    """
    𝜀̃ₛ::Float64
    𝜀̃∞::Float64
    τ::Float64
    function DebyeMaterial(𝜀̃ₛ::Real, 𝜀̃∞::Real, τ::Real)
        𝜀̃ₛ > 0  || error("Debye static permittivity 𝜀̃ₛ must be positive, got $𝜀̃ₛ")
        𝜀̃∞ > 0  || error("Debye high-frequency permittivity 𝜀̃∞ must be positive, got $𝜀̃∞")
        τ > 0   || error("Debye relaxation time τ must be positive, got $τ")
        new(Float64(𝜀̃ₛ), Float64(𝜀̃∞), Float64(τ))
    end
end

# The permittivity of a DebyeMaterial is computed using the formula given in the docstring, which accounts for the static and high-frequency permittivities and the relaxation time.
function 𝜀̃(material::DebyeMaterial, f::Real)
    ω = 2π * f
    δ𝜀̃ = material.𝜀̃ₛ - material.𝜀̃∞
    material.𝜀̃∞ + δ𝜀̃ / (1.0 + 1im*ω*material.τ)
end


# ---- SampledDataMaterial Definition ----
struct SampledDataMaterial <: AbstractMaterial
    """
    SampledDataMaterial(wavelengths::AbstractVector{<:Real},
                        n::AbstractVector{<:Real},
                        k::AbstractVector{<:Real};
                        interp_kind=Gridded(Linear()))
    Material defined by sampled refractive index and extinction coefficient data as a function of wavelength.
    - wavelengths: Vector of wavelengths in nanometres (must be positive and in ascending order)
    - n: Vector of refractive index values corresponding to the wavelengths (must be non-negative)
    - k: Vector of extinction coefficient values corresponding to the wavelengths (must be non-negative)
    - interp_kind: Interpolation method for n and k (default is linear on the non-uniform
                   wavelength grid; pass e.g. `Gridded(Linear())` or another `Gridded` scheme)
    The permittivity is computed by interpolating n and k at the wavelength corresponding to the given frequency f, and then using 𝜀̃ = (n + i*k)^2.
    Note: 
        * The wavelengths must be in ascending order and correspond to the n and k values.
        * For frequencies outside the range of the provided wavelengths, the interpolation will extrapolate based on the chosen method.
    """
    wavelengths::Vector{Float64}
    n_data::Vector{Float64}
    k_data::Vector{Float64}
    interp_n
    interp_k

    function SampledDataMaterial(wavelengths::AbstractVector{<:Real},
                                n::AbstractVector{<:Real},
                                k::AbstractVector{<:Real};
                                interp_kind=Gridded(Linear()))
        (length(wavelengths) == length(n) == length(k)) ||
            error("wavelengths, n, k must have the same length")
        wl = Float64.(wavelengths)
        issorted(wl) || error("wavelengths must be in ascending order")
        n_vals = Float64.(n)
        k_vals = Float64.(k)
        itp_n = interpolate((wl,), n_vals, interp_kind)
        itp_k = interpolate((wl,), k_vals, interp_kind)
        new(wl, n_vals, k_vals, itp_n, itp_k)
    end

    function SampledDataMaterial(datafile::String; interp_kind=Gridded(Linear()))
        # Load data from file
        data = readdlm(datafile)
        if size(data, 2) != 3
            error("Data file must have exactly three columns: wavelength, n, k")
        end
        wavelengths = data[:, 1]
        n = data[:, 2]
        k = data[:, 3]
        SampledDataMaterial(wavelengths, n, k; interp_kind=interp_kind)
    end

end

# The permittivity of a SampledDataMaterial is computed by first converting the frequency f to wavelength, then interpolating n and k at that wavelength, and finally using the formula 𝜀̃ = (n + i*k)^2.
function 𝜀̃(material::SampledDataMaterial, f::Real)
    λ_m = 299792458.0 / f
    λ_nm = λ_m * 1e9
    n = material.interp_n(λ_nm)
    k = material.interp_k(λ_nm)
    Complex{Float64}(n + 1im*k)^2
end

function ñ(material::SampledDataMaterial, f::Real)
    λ_m = 299792458.0 / f
    λ_nm = λ_m * 1e9
    n = material.interp_n(λ_nm)
    k = material.interp_k(λ_nm)
    Complex{Float64}(n, k)
end