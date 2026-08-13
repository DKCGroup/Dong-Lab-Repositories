using LinearAlgebra

# CMT — tight-binding coupled-mode theory for twisted/layered photonic stacks. Each resonant
# mode j has a complex on-site frequency ω₀ⱼ + iκ₀ⱼ (resonance + decay). Modes couple through a
# characteristic ("Hamiltonian") matrix H whose off-diagonal entries are the intralayer and
# interlayer coupling coefficients (Eqs. (S3)-(S7), (S28)-(S33) of the supplement). The complex
# eigenvalues of H are the complex eigenfrequencies of the system: Re gives the resonance
# frequency, Im the linewidth/decay, and Q = Re / (2·|Im|).
#
# For periodic systems the couplings carry Bloch phases Σ_δ exp(-i k·δ) (Eqs. (S12)-(S16)), so
# H = H(k) and sweeping k yields the photonic band structure. Parameters (ω₀, κ₀, coupling g)
# can be supplied directly or extracted from first-principles RCWA simulations — the "fusion"
# of the two methods — using `extract_resonance` and `extract_coupling`.
if !isdefined(@__MODULE__, :c₀)
    const c₀ = 299792458.0
end


# ---- Coupled-mode characteristic matrix ----
# Assemble H from on-site complex frequencies (placed on the diagonal) and a coupling matrix
# (off-diagonal). `onsite[j]` should be ω₀ⱼ + iκ₀ⱼ; `coupling[j,k]` is the coupling from mode k
# into mode j (its diagonal is ignored).
function cmt_matrix(onsite::AbstractVector, coupling::AbstractMatrix)
    n = length(onsite)
    size(coupling) == (n, n) || error("coupling must be $n×$n to match onsite")
    H = Matrix{Complex{Float64}}(coupling)
    for j in 1:n
        H[j, j] = onsite[j]
    end
    H
end

# Complex eigenfrequencies of a characteristic matrix, sorted by ascending resonance (real
# part). Each eigenvalue Ω = ω + iγ encodes a resonance at ω with linewidth/decay γ.
function eigenfrequencies(H::AbstractMatrix)
    Ω = eigvals(Matrix{Complex{Float64}}(H))
    sort!(Ω; by=real)
    Ω
end

# Resonance frequencies (real parts) of the complex eigenfrequencies.
resonances(Ω::AbstractVector) = real.(Ω)

# Linewidths / decay rates (absolute imaginary parts) of the complex eigenfrequencies.
linewidths(Ω::AbstractVector) = abs.(imag.(Ω))

# Quality factors Q = ω / (2·|γ|) for each complex eigenfrequency.
function qfactors(Ω::AbstractVector)
    [imag(ω) == 0 ? Inf : abs(real(ω)) / (2 * abs(imag(ω))) for ω in Ω]
end


# ---- Bloch (k-dependent) coupled-mode matrix ----
struct CMTCoupling
    """
    CMTCoupling(i, j, g, δ)
    A directed coupling that adds g·exp(-i k·δ) to H[i, j] at Bloch wavevector k.
    - i, j: Mode indices coupled (into i from j)
    - g: Coupling coefficient (complex)
    - δ: Real-space displacement (δx, δy) carried by this hopping term
    """
    i::Int
    j::Int
    g::Complex{Float64}
    δ::NTuple{2,Float64}
end

# Build the Bloch characteristic matrix H(k) from on-site complex frequencies and a list of
# directed couplings with displacements (Eqs. (S12)-(S14)).
function bloch_matrix(onsite::AbstractVector, couplings::AbstractVector{CMTCoupling},
                      k::NTuple{2,<:Real})
    n = length(onsite)
    H = zeros(Complex{Float64}, n, n)
    for j in 1:n
        H[j, j] = onsite[j]
    end
    for c in couplings
        (1 <= c.i <= n && 1 <= c.j <= n) || error("coupling index out of range")
        H[c.i, c.j] += c.g * cis(-(k[1] * c.δ[1] + k[2] * c.δ[2]))
    end
    H
end

# Photonic band structure: complex eigenfrequencies at each Bloch wavevector in `ks`.
# Returns a matrix of size (length(ks), n) with eigenfrequencies sorted per k.
function bands(onsite::AbstractVector, couplings::AbstractVector{CMTCoupling},
               ks::AbstractVector{<:NTuple{2,<:Real}})
    n = length(onsite)
    out = Matrix{Complex{Float64}}(undef, length(ks), n)
    for (p, k) in enumerate(ks)
        out[p, :] = eigenfrequencies(bloch_matrix(onsite, couplings, k))
    end
    out
end


# ---- Fusion: extract CMT parameters from RCWA spectra ----
# Locate a single isolated resonance in a spectrum (e.g. transmittance or reflectance vs.
# frequency) and return its centre frequency f₀ and half-width-at-half-maximum κ (the CMT
# decay rate). Set `dip=true` for a transmission dip (resonant absorption/reflection) or
# `dip=false` for a peak. Frequencies `fs` must be sorted ascending.
function extract_resonance(fs::AbstractVector{<:Real}, spectrum::AbstractVector{<:Real};
                           dip::Bool=true)
    length(fs) == length(spectrum) || error("fs and spectrum must have equal length")
    length(fs) >= 3 || error("need at least 3 sample points")
    s = dip ? -collect(float.(spectrum)) : collect(float.(spectrum))   # turn dip into a peak
    ipk = argmax(s)
    f0 = fs[ipk]
    base = (s[1] + s[end]) / 2          # background level away from resonance
    half = base + (s[ipk] - base) / 2   # half-maximum level (above background)
    # walk left and right to the half-maximum crossings, linearly interpolating
    fl = fs[1]
    for i in ipk:-1:2
        if s[i] >= half >= s[i-1]
            fl = fs[i-1] + (fs[i] - fs[i-1]) * (half - s[i-1]) / (s[i] - s[i-1])
            break
        end
    end
    fr = fs[end]
    for i in ipk:(length(fs)-1)
        if s[i] >= half >= s[i+1]
            fr = fs[i] + (fs[i+1] - fs[i]) * (s[i] - half) / (s[i] - s[i+1])
            break
        end
    end
    κ = (fr - fl) / 2
    (f0=f0, κ=κ)
end

# Extract a coupling coefficient from the frequency splitting of two hybridised modes (e.g. the
# symmetric/antisymmetric pair of an RCWA-computed coupled bilayer): g = |f₊ − f₋| / 2 (Eq. (S5)).
extract_coupling(f_plus::Real, f_minus::Real) = abs(f_plus - f_minus) / 2
