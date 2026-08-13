# api.jl — cross-language entry points for the IR-metasurface workflow.
#
# These functions are the Julia side of the Python API (see python/gestalt_metasurface.py).
# Design rules for smooth juliacall interop:
#   • accept/return only plain data (Vector/Matrix/Float64/tuple/NamedTuple) — no Plots;
#   • grayscale images are passed as 2-D real matrices (values 0–255), any numeric element
#     type is accepted and converted internally;
#   • the VO2 phase is selected by a string ("insulating" / "metallic");
#   • all keyword-argument names are ASCII so Python can forward them directly.
#
# Include after `metasurface_ir.jl` (or after `using GESTALT`).

if !isdefined(@__MODULE__, :vo2_insulating_sampled)
    include(joinpath(@__DIR__, "metasurface_ir.jl"))
end

# Resolve a phase string to a VO2 material. `sampled=true` uses the measured n,k table;
# `sampled=false` falls back to the Lorentz (insulating) / Drude (metallic) models.
function _phase_material(phase::String; sampled::Bool=true)
    p = lowercase(strip(phase))
    if p in ("insulating", "ins", "i", "dielectric", "monoclinic")
        return sampled ? vo2_insulating_sampled() : vo2_insulating_material()
    elseif p in ("metallic", "met", "m", "metal", "rutile")
        return sampled ? vo2_metallic_sampled() : vo2_metallic_material()
    else
        error("unknown VO2 phase '$phase'; use \"insulating\" or \"metallic\"")
    end
end

# VO2 complex refractive index ñ = n + i·k over the given wavelengths (metres).
# Returns (wavelength, n, k) — handy for verifying/plotting the material dispersion.
function vo2_dispersion(phase::String; lambdas=range(8e-6, 14e-6; length=200),
                        sampled::Bool=true)
    mat = _phase_material(phase; sampled=sampled)
    λ = collect(Float64.(lambdas))
    nk = [ñ(mat, c₀ / λi) for λi in λ]
    (wavelength=λ, n=real.(nk), k=imag.(nk))
end

# Core routine: grayscale image → R/T/A spectrum.
#   img  : N×N real matrix, values 0–255 (0 = air, v > t_min_gray = VO2 of height v/255·max_t)
#   phase: "insulating" or "metallic"
# Returns (wavelength, R, T, A) each a Vector{Float64} of length n_lambda, wavelength in metres.
function metasurface_spectrum(img::AbstractMatrix, phase::String;
                              p::Real=5e-6, max_t::Real=4e-6, t_min_gray::Real=50,
                              N::Integer=6, lambda_min::Real=8e-6, lambda_max::Real=14e-6,
                              n_lambda::Integer=200, theta::Real=0.0, sampled::Bool=true,
                              d_al::Real=200e-9)
    n_lambda >= 2 || error("n_lambda must be >= 2, got $n_lambda")
    imgf = Float64.(img)
    λs = collect(range(lambda_min, lambda_max; length=n_lambda))
    mat = _phase_material(phase; sampled=sampled)
    al = aluminum_material()
    grid, layers = build_stack(imgf, mat, al; p=p, max_t=max_t, t_min_gray=t_min_gray,
                               N=N, d_al=d_al)
    R, T, A = sweep_spectrum(layers, grid, λs; θ=theta)
    (wavelength=λs, R=R, T=T, A=A)
end

# Convenience: shape name + phase → R/T/A spectrum in one call.
function metasurface_shape_spectrum(shape::String, phase::String; kwargs...)
    metasurface_spectrum(shape_image(shape), phase; kwargs...)
end
