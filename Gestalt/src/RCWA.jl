using LinearAlgebra

# Layer.jl is included (not `using Layer`) to match the flat layout of this repo; it in turn
# loads Material.jl. Load it before RCWA references AbstractLayer, material, thickness, 𝜀̃,
# and the periodic-pattern accessors (g_vectors, ε̂).
if !isdefined(@__MODULE__, :AbstractLayer)
    include(joinpath(@__DIR__, "Layer.jl"))
end

# Speed of light in vacuum (m/s), matching the constant used in Material.jl / TMM.jl.
if !isdefined(@__MODULE__, :c₀)
    const c₀ = 299792458.0
end


# ---- Reciprocal-space grid (plane-wave truncation) ----
struct RCWAGrid
    """
    RCWAGrid(b1, b2; Nx::Integer=0, Ny::Integer=0)
    Truncated reciprocal-lattice basis used to expand all fields in plane waves
    (Eqs. (1)-(3) of the supplement): g = m·b1 + n·b2 with m ∈ -Nx..Nx, n ∈ -Ny..Ny.
    - b1, b2: Reciprocal lattice vectors (gx, gy) of the layer lattice (rad / length)
    - Nx, Ny: Truncation orders along b1 and b2 (the total harmonic count is (2Nx+1)(2Ny+1))
    The harmonics are stored together with their integer orders (m, n) so that convolution
    matrices [εr]_{p,q} = ε̂(g_p − g_q) can be looked up exactly by order differences.
    """
    b1::NTuple{2,Float64}
    b2::NTuple{2,Float64}
    orders::Vector{NTuple{2,Int}}
    gx::Vector{Float64}
    gy::Vector{Float64}
    center::Int
    function RCWAGrid(b1::NTuple{2,<:Real}, b2::NTuple{2,<:Real}; Nx::Integer=0, Ny::Integer=0)
        Nx >= 0 || error("Nx must be non-negative, got $Nx")
        Ny >= 0 || error("Ny must be non-negative, got $Ny")
        b1f = NTuple{2,Float64}(b1)
        b2f = NTuple{2,Float64}(b2)
        orders = NTuple{2,Int}[]
        gx = Float64[]
        gy = Float64[]
        center = 0
        for n in (-Ny):Ny, m in (-Nx):Nx
            push!(orders, (m, n))
            push!(gx, m * b1f[1] + n * b2f[1])
            push!(gy, m * b1f[2] + n * b2f[2])
            (m == 0 && n == 0) && (center = length(orders))
        end
        center != 0 || error("zero order (0, 0) missing from grid")
        new(b1f, b2f, orders, gx, gy, center)
    end
end

# Number of retained harmonics (plane-wave bases) in the grid.
nharmonics(grid::RCWAGrid) = length(grid.orders)


# ---- Convolution matrices [εr] ----
# Convolution (Toeplitz-like) matrix of the relative permittivity on the grid:
# [εr]_{p,q} = ε̂(g_p − g_q). For a homogeneous layer this reduces to ε·I; for a periodic
# layer the Fourier amplitudes are taken from its FourierPeriodicPattern, keyed by integer
# order differences (m_p − m_q, n_p − n_q).
function convolution_matrix(layer::AbstractLayer, grid::RCWAGrid, f::Real)
    NH = nharmonics(grid)
    if layer isa HomogeneousLayer
        return Matrix{Complex{Float64}}(𝜀̃(material(layer), f) * I, NH, NH)
    elseif layer isa PeriodicLayer && layer.pattern isa FourierPeriodicPattern
        coeffs = _fourier_order_map(layer.pattern, grid, f)
        ER = zeros(Complex{Float64}, NH, NH)
        for p in 1:NH, q in 1:NH
            Δm = grid.orders[p][1] - grid.orders[q][1]
            Δn = grid.orders[p][2] - grid.orders[q][2]
            ER[p, q] = get(coeffs, (Δm, Δn), 0.0 + 0.0im)
        end
        return ER
    else
        error("convolution_matrix not supported for $(typeof(layer)) with this pattern; " *
              "use HomogeneousLayer or PeriodicLayer{FourierPeriodicPattern}")
    end
end

# Map a FourierPeriodicPattern's Fourier amplitudes onto integer grid orders (m, n) by
# matching each stored g-vector to the closest m·b1 + n·b2. Coefficients given as functions
# of f are evaluated at f.
function _fourier_order_map(pattern::FourierPeriodicPattern, grid::RCWAGrid, f::Real)
    coeffs = Dict{NTuple{2,Int},Complex{Float64}}()
    B = [grid.b1[1] grid.b2[1]; grid.b1[2] grid.b2[2]]
    for (g, c) in zip(pattern.g_vectors, pattern.ε̂_coefficients)
        mn = B \ [g[1]; g[2]]
        m = round(Int, mn[1])
        n = round(Int, mn[2])
        (abs(mn[1] - m) < 1e-6 && abs(mn[2] - n) < 1e-6) ||
            error("Fourier g-vector $g is not an integer combination of (b1, b2)")
        coeffs[(m, n)] = _eval_ε̂_coeff(c, f)
    end
    coeffs
end


# ---- Normalised in-plane wavevectors ----
# Diagonal (per-harmonic) normalised in-plane wavevectors K̃x, K̃y = (kinc + g) / k0, where
# the incidence medium fixes kinc through the angles (θ, φ). Returned as plain vectors.
function _normalized_k(grid::RCWAGrid, k0::Real, n_inc, θ::Real, φ::Real)
    kx0 = real(n_inc) * sin(θ) * cos(φ)
    ky0 = real(n_inc) * sin(θ) * sin(φ)
    kx = [kx0 + grid.gx[i] / k0 for i in 1:nharmonics(grid)]
    ky = [ky0 + grid.gy[i] / k0 for i in 1:nharmonics(grid)]
    kx, ky
end


# ---- Block helpers ----
# Assemble a 2×2 block matrix [A B; C D] as a dense complex matrix.
function _block(A, B, C, D)
    [Matrix{Complex{Float64}}(A) Matrix{Complex{Float64}}(B);
     Matrix{Complex{Float64}}(C) Matrix{Complex{Float64}}(D)]
end


# Principal-like square root with a fixed branch (Im ≥ 0; tie-break Re ≥ 0). This selects the
# physical out-of-plane wavevector kz of a forward wave that decays along +z under the e^{-iωt}
# time convention (matching Material.jl / TMM.jl), so passive media never act as gain media.
function _sqrt_kz(z)
    s = sqrt(Complex{Float64}(z))
    (imag(s) < 0 || (imag(s) == 0 && real(s) < 0)) ? -s : s
end

# ---- Homogeneous-region eigenmodes ----
# Out-of-plane wavevectors KZ of a homogeneous region of permittivity ε (μ = 1), on the
# Im(kz) ≥ 0 branch so that evanescent orders decay along +z (e^{-iωt} convention).
function _region_kz(ε, kx::AbstractVector, ky::AbstractVector)
    [_sqrt_kz(ε - kx[i]^2 - ky[i]^2) for i in eachindex(kx)]
end

# Eigenmodes (W, V) of a homogeneous region with permittivity ε (μ = 1). W = I and the
# magnetic modes follow from V = Q·λ⁻¹ with the analytic region Q-matrix, where the modal
# eigenvalue is λ = -i·kz. This matches the layer-mode branch (Re λ ≥ 0; propagating modes
# have Im λ ≤ 0), which is essential for correct energy accounting in lossy stacks.
function _region_modes(ε, kx::AbstractVector, ky::AbstractVector)
    NH = length(kx)
    Inh = Diagonal(ones(Complex{Float64}, NH))
    KX = Diagonal(Complex{Float64}.(kx))
    KY = Diagonal(Complex{Float64}.(ky))
    kz = _region_kz(ε, kx, ky)
    Q = _block(KX * KY, ε * Inh - KX * KX, KY * KY - ε * Inh, -KX * KY)
    W = Matrix{Complex{Float64}}(I, 2NH, 2NH)
    λ = -1im .* vcat(kz, kz)
    # Floor tiny |λ| (Wood anomalies / grazing: some order has kz ≈ 0) so 1/λ stays finite.
    λsafe = [abs(λi) < 1e-12 ? Complex{Float64}(1e-12, 0) : λi for λi in λ]
    V = Q * Diagonal(1.0 ./ λsafe)
    W, V, kz
end


# ---- Layer eigenmodes ----
# Eigenmodes of a (possibly patterned) layer from Ω² = P·Q (Eqs. (14)-(16), μ = 1). Returns
# the electric-field modes W, the magnetic-field modes V = Q·W·λ⁻¹, and the per-mode λ
# (chosen with Re λ ≥ 0 so that exp(−λ·k0·d) decays through the layer).
function _layer_modes(ER::AbstractMatrix, kx::AbstractVector, ky::AbstractVector)
    NH = length(kx)
    Inh = Diagonal(ones(Complex{Float64}, NH))
    KX = Diagonal(Complex{Float64}.(kx))
    KY = Diagonal(Complex{Float64}.(ky))
    ERi = inv(ER)
    P = _block(KX * ERi * KY, Inh - KX * ERi * KX, KY * ERi * KY - Inh, -KY * ERi * KX)
    Q = _block(KX * KY, ER - KX * KX, KY * KY - ER, -KY * KX)
    F = eigen(P * Q)
    λ = sqrt.(Complex{Float64}.(F.values))
    for i in eachindex(λ)
        # Forward-decaying branch: Re λ ≥ 0 (no growth through the layer); for purely
        # propagating modes (Re λ = 0) take Im λ ≤ 0 to match the region/gap convention.
        (real(λ[i]) < 0 || (real(λ[i]) == 0 && imag(λ[i]) > 0)) && (λ[i] = -λ[i])
    end
    W = F.vectors
    V = Q * W * Diagonal(1.0 ./ λ)
    W, V, λ
end

# Analytic eigenmodes of a laterally homogeneous layer of permittivity ε: identical to the
# region modes (W = I, V = Q·λ⁻¹, λ = -i·kz on the Im kz ≥ 0 branch). Avoids the dense
# eigendecomposition of `_layer_modes` for unpatterned inner layers (spacers, metal films),
# giving the same layer scattering matrix at a fraction of the cost.
function _homogeneous_modes(ε, kx::AbstractVector, ky::AbstractVector)
    W, V, kz = _region_modes(ε, kx, ky)
    λ = -1im .* vcat(kz, kz)
    W, V, λ
end


# ---- Scattering matrices ----
struct SMatrix
    """
    SMatrix(S11, S12, S21, S22)
    Block scattering matrix relating incoming to outgoing mode amplitudes on both sides,
    b = S·a, used with the Redheffer star product to cascade interfaces and layers.
    """
    S11::Matrix{Complex{Float64}}
    S12::Matrix{Complex{Float64}}
    S21::Matrix{Complex{Float64}}
    S22::Matrix{Complex{Float64}}
end

# Redheffer star product S = SA ⋆ SB cascading two scattering matrices.
function star(SA::SMatrix, SB::SMatrix)
    n = size(SA.S11, 1)
    Id = Matrix{Complex{Float64}}(I, n, n)
    D = SA.S12 * inv(Id - SB.S11 * SA.S22)
    F = SB.S21 * inv(Id - SA.S22 * SB.S11)
    S11 = SA.S11 + D * SB.S11 * SA.S21
    S12 = D * SB.S12
    S21 = F * SA.S21
    S22 = SB.S22 + F * SA.S22 * SB.S12
    SMatrix(S11, S12, S21, S22)
end

# Scattering matrix of a single layer referenced to the gap medium (Wg, Vg), with modal
# propagation X = exp(−λ·k0·d).
function _layer_smatrix(W, V, λ, Wg, Vg, k0::Real, d::Real)
    A = W \ Wg + V \ Vg
    B = W \ Wg - V \ Vg
    Ai = inv(A)
    X = Diagonal(exp.(-λ .* (k0 * d)))
    XB = X * B
    factor = inv(A - XB * Ai * XB)
    S11 = factor * (XB * Ai * X * A - B)
    S12 = factor * X * (A - B * Ai * B)
    SMatrix(S11, S12, S12, S11)
end

# Scattering matrix connecting the external (semi-infinite) reflection region to the gap.
function _reflection_smatrix(Wref, Vref, Wg, Vg)
    A = Wg \ Wref + Vg \ Vref
    B = Wg \ Wref - Vg \ Vref
    Ai = inv(A)
    S11 = -Ai * B
    S12 = 2 * Ai
    S21 = 0.5 * (A - B * Ai * B)
    S22 = B * Ai
    SMatrix(S11, S12, S21, S22)
end

# Scattering matrix connecting the gap to the external (semi-infinite) transmission region.
function _transmission_smatrix(Wtrn, Vtrn, Wg, Vg)
    A = Wg \ Wtrn + Vg \ Vtrn
    B = Wg \ Wtrn - Vg \ Vtrn
    Ai = inv(A)
    S11 = B * Ai
    S12 = 0.5 * (A - B * Ai * B)
    S21 = 2 * Ai
    S22 = -Ai * B
    SMatrix(S11, S12, S21, S22)
end


# ---- 2×2 block-diagonal scattering matrices (per-harmonic fast path) ----------------------
# For a laterally homogeneous region (including the ε = 1 gap medium) the eigenmode basis is
# block-diagonal: W = I and V = Q·λ⁻¹ with Q built from diagonal Kx, Ky, so V splits into an
# independent 2×2 block per harmonic (the (x, y) polarization subspace). Consequently every
# homogeneous layer and every homogeneous↔homogeneous interface has a 2×2 block-diagonal
# scattering matrix, and the Redheffer star product of two such matrices is again 2×2
# block-diagonal. Cascading M homogeneous layers therefore costs O(M·NH) 2×2 operations
# instead of O(M·NH³) dense ones; only patterned layers break the decoupling and require the
# dense path below.

struct BlockDiagS
    S11::Array{Complex{Float64},3}
    S12::Array{Complex{Float64},3}
    S21::Array{Complex{Float64},3}
    S22::Array{Complex{Float64},3}
    NH::Int
end

const _I2 = Matrix{Complex{Float64}}(I, 2, 2)

# (V, kz) of a homogeneous region of permittivity ε on the per-harmonic 2×2 basis.
# V = Q·λ⁻¹ uses the same λ-flooring as `_region_modes`; kz is the raw branch-selected
# out-of-plane wavevector (Im kz ≥ 0).
function _region_blockdiag(ε, kx::AbstractVector, ky::AbstractVector)
    NH = length(kx)
    V = Array{Complex{Float64}}(undef, 2, 2, NH)
    kz = Vector{Complex{Float64}}(undef, NH)
    for i in 1:NH
        kzi = _sqrt_kz(ε - kx[i]^2 - ky[i]^2)
        kz[i] = kzi
        λi = -1im * kzi
        λsafe = abs(λi) < 1e-12 ? Complex{Float64}(1e-12, 0) : λi
        invλ = 1.0 / λsafe
        q11 = kx[i] * ky[i]; q12 = ε - kx[i]^2
        q21 = ky[i]^2 - ε; q22 = -kx[i] * ky[i]
        V[1, 1, i] = q11 * invλ; V[1, 2, i] = q12 * invλ
        V[2, 1, i] = q21 * invλ; V[2, 2, i] = q22 * invλ
    end
    V, kz
end

# Redheffer star product of two 2×2 block-diagonal scattering matrices (per harmonic).
function star(SA::BlockDiagS, SB::BlockDiagS)
    NH = SA.NH
    NH == SB.NH || error("BlockDiagS size mismatch")
    S11 = Array{Complex{Float64}}(undef, 2, 2, NH)
    S12 = similar(S11); S21 = similar(S11); S22 = similar(S11)
    @inbounds for i in 1:NH
        a11 = SA.S11[:, :, i]; a12 = SA.S12[:, :, i]
        a21 = SA.S21[:, :, i]; a22 = SA.S22[:, :, i]
        b11 = SB.S11[:, :, i]; b12 = SB.S12[:, :, i]
        b21 = SB.S21[:, :, i]; b22 = SB.S22[:, :, i]
        D = a12 * inv(_I2 - b11 * a22)
        F = b21 * inv(_I2 - a22 * b11)
        S11[:, :, i] = a11 + D * b11 * a21
        S12[:, :, i] = D * b12
        S21[:, :, i] = F * a21
        S22[:, :, i] = b22 + F * a22 * b12
    end
    BlockDiagS(S11, S12, S21, S22, NH)
end

# Dense ⇄ block-diagonal conversions.
function _bd_to_dense(B::BlockDiagS)
    NH = B.NH
    N = 2NH
    S11 = zeros(Complex{Float64}, N, N)
    S12 = zeros(Complex{Float64}, N, N)
    S21 = zeros(Complex{Float64}, N, N)
    S22 = zeros(Complex{Float64}, N, N)
    for i in 1:NH
        i2 = NH + i
        S11[i, i] = B.S11[1, 1, i]; S11[i, i2] = B.S11[1, 2, i]
        S11[i2, i] = B.S11[2, 1, i]; S11[i2, i2] = B.S11[2, 2, i]
        S12[i, i] = B.S12[1, 1, i]; S12[i, i2] = B.S12[1, 2, i]
        S12[i2, i] = B.S12[2, 1, i]; S12[i2, i2] = B.S12[2, 2, i]
        S21[i, i] = B.S21[1, 1, i]; S21[i, i2] = B.S21[1, 2, i]
        S21[i2, i] = B.S21[2, 1, i]; S21[i2, i2] = B.S21[2, 2, i]
        S22[i, i] = B.S22[1, 1, i]; S22[i, i2] = B.S22[1, 2, i]
        S22[i2, i] = B.S22[2, 1, i]; S22[i2, i2] = B.S22[2, 2, i]
    end
    SMatrix(S11, S12, S21, S22)
end

# Mixed star products: promote to dense only at the (patterned) dense boundary.
star(SA::BlockDiagS, SB::SMatrix) = star(_bd_to_dense(SA), SB)
star(SA::SMatrix, SB::BlockDiagS) = star(SA, _bd_to_dense(SB))

# 2×2 block-diagonal scattering matrix of a homogeneous layer (permittivity ε, thickness d)
# referenced to the gap medium (Vg block). Vl and kzl come from `_region_blockdiag`.
function _layer_smatrix_blockdiag(Vl::Array{Complex{Float64},3}, kzl::AbstractVector,
                                  k0::Real, d::Real, Vg::Array{Complex{Float64},3})
    NH = length(kzl)
    S11 = Array{Complex{Float64}}(undef, 2, 2, NH)
    S12 = similar(S11); S21 = similar(S11); S22 = similar(S11)
    for i in 1:NH
        Vli = Vl[:, :, i]; Vgi = Vg[:, :, i]
        # A = W\Wg + V\Vg = I + V_layer⁻¹·V_gap,  B = I − V_layer⁻¹·V_gap
        A = _I2 + Vli \ Vgi
        B = _I2 - Vli \ Vgi
        Ai = inv(A)
        # X = e·I with e = exp(−λ·k0·d) = exp(i·kz·k0·d)  (x and y modes share λ = −i·kz)
        e = exp(1im * kzl[i] * (k0 * d))
        X = e * _I2
        XB = X * B
        factor = inv(A - XB * Ai * XB)
        S11[:, :, i] = factor * (XB * Ai * X * A - B)
        S12[:, :, i] = factor * X * (A - B * Ai * B)
        S21[:, :, i] = S12[:, :, i]
        S22[:, :, i] = S11[:, :, i]
    end
    BlockDiagS(S11, S12, S21, S22, NH)
end

# 2×2 block-diagonal reflection (ref half-space → gap) scattering matrix.
function _reflection_smatrix_blockdiag(Vref::Array{Complex{Float64},3},
                                       Vg::Array{Complex{Float64},3})
    NH = size(Vref, 3)
    S11 = Array{Complex{Float64}}(undef, 2, 2, NH)
    S12 = similar(S11); S21 = similar(S11); S22 = similar(S11)
    for i in 1:NH
        A = _I2 + Vg[:, :, i] \ Vref[:, :, i]    # Wg\Wref = I
        B = _I2 - Vg[:, :, i] \ Vref[:, :, i]
        Ai = inv(A)
        S11[:, :, i] = -Ai * B
        S12[:, :, i] = 2 * Ai
        S21[:, :, i] = 0.5 * (A - B * Ai * B)
        S22[:, :, i] = B * Ai
    end
    BlockDiagS(S11, S12, S21, S22, NH)
end

# 2×2 block-diagonal transmission (gap → trn half-space) scattering matrix.
function _transmission_smatrix_blockdiag(Vtrn::Array{Complex{Float64},3},
                                         Vg::Array{Complex{Float64},3})
    NH = size(Vtrn, 3)
    S11 = Array{Complex{Float64}}(undef, 2, 2, NH)
    S12 = similar(S11); S21 = similar(S11); S22 = similar(S11)
    for i in 1:NH
        A = _I2 + Vg[:, :, i] \ Vtrn[:, :, i]
        B = _I2 - Vg[:, :, i] \ Vtrn[:, :, i]
        Ai = inv(A)
        S11[:, :, i] = B * Ai
        S12[:, :, i] = 0.5 * (A - B * Ai * B)
        S21[:, :, i] = 2 * Ai
        S22[:, :, i] = -Ai * B
    end
    BlockDiagS(S11, S12, S21, S22, NH)
end


# ---- RCWAResult Definition ----
struct RCWAResult
    """
    RCWAResult(orders, R, T, Rtot, Ttot, A)
    Diffraction efficiencies of an RCWA calculation for a single (f, θ, φ, polarization).
    - orders: Integer diffraction orders (m, n) of the grid, in the same order as R and T
    - R, T: Per-order reflectance and transmittance (real, ≥ 0 for propagating orders)
    - Rtot, Ttot: Summed reflectance and transmittance over all orders
    - A: Absorptance, 1 − Rtot − Ttot
    """
    orders::Vector{NTuple{2,Int}}
    R::Vector{Float64}
    T::Vector{Float64}
    Rtot::Float64
    Ttot::Float64
    A::Float64
end


# ---- Scattering-matrix assembly at fixed in-plane Bloch momentum ----
# Absolute in-plane wavevector (kx_abs, ky_abs) [rad/m] → normalised harmonics
# kx = (kx_abs + g_x)/k0, ky = (ky_abs + g_y)/k0. Works above and below the light line
# (below: |k∥| > n k0 ⇒ all exterior orders are evanescent).
function _kxky_from_kabs(grid, k0::Real, kx_abs::Real, ky_abs::Real)
    NH = nharmonics(grid)
    kx = [kx_abs / k0 + grid.gx[i] / k0 for i in 1:NH]
    ky = [ky_abs / k0 + grid.gy[i] / k0 for i in 1:NH]
    kx, ky
end

# Cascade the full S-matrix of a layer stack at frequency f and Bloch momentum (kx_abs, ky_abs).
# Returns (SG, kx, ky, kz_ref, kz_trn, k0, n_ref). Used by `solve` and by guided-mode maps.
# Homogeneous layers and half-space interfaces are cascaded on the 2×2 block-diagonal fast
# path (O(NH)); only patterned layers incur dense O(NH³) work.
function smatrix_k(layers::AbstractVector{<:AbstractLayer}, grid::RCWAGrid, f::Real,
                   kx_abs::Real, ky_abs::Real; deriv::Symbol=:spectral)
    N = length(layers)
    N >= 2 || error("need at least incidence and transmission half-spaces, got $N")
    (layers[1] isa HomogeneousLayer && layers[N] isa HomogeneousLayer) ||
        error("first and last layers must be HomogeneousLayer (semi-infinite regions)")

    k0 = 2π * f / c₀
    ε_ref = 𝜀̃(material(layers[1]), f)
    ε_trn = 𝜀̃(material(layers[N]), f)
    n_ref = sqrt(ε_ref)
    kx, ky = _kxky_from_kabs(grid, k0, kx_abs, ky_abs)

    # Dense gap modes (needed only by patterned layers).
    Wg, Vg, _ = _region_modes(1.0, kx, ky)
    # Block-diagonal region data (used by all homogeneous layers/interfaces).
    Vg_bd, _ = _region_blockdiag(1.0, kx, ky)
    Vref_bd, kz_ref = _region_blockdiag(ε_ref, kx, ky)
    Vtrn_bd, kz_trn = _region_blockdiag(ε_trn, kx, ky)

    SG = _reflection_smatrix_blockdiag(Vref_bd, Vg_bd)   # BlockDiagS
    for i in 2:(N - 1)
        layer = layers[i]
        d = thickness(layer)
        isfinite(d) || error("inner layer $i must have finite thickness, got $d")
        if layer isa HomogeneousLayer
            Vl_bd, kzl = _region_blockdiag(𝜀̃(material(layer), f), kx, ky)
            SG = star(SG, _layer_smatrix_blockdiag(Vl_bd, kzl, k0, d, Vg_bd))
        elseif layer isa PeriodicLayer && pattern(layer) isa PixelPeriodicPattern
            pat = pattern(layer)
            ε_in = 𝜀̃(pat.mat_in, f)
            ε_out = 𝜀̃(pat.mat_out, f)
            Sl = layer_smatrix_realspace(pat.indicator, ε_in, ε_out, grid,
                                         k0, kx_abs, ky_abs, d, pat.period,
                                         size(pat.indicator, 1); deriv=deriv)
            SG = star(SG, Sl)
        else
            ER = convolution_matrix(layer, grid, f)
            W, V, λ = _layer_modes(ER, kx, ky)
            SG = star(SG, _layer_smatrix(W, V, λ, Wg, Vg, k0, d))
        end
    end
    SG = star(SG, _transmission_smatrix_blockdiag(Vtrn_bd, Vg_bd))
    SG isa BlockDiagS && (SG = _bd_to_dense(SG))
    (SG, kx, ky, kz_ref, kz_trn, k0, n_ref)
end

# Resonance indicator from the scattering matrix: large near guided / quasi-guided poles.
# Uses log₁₀ of the spectral norm of S₁₁ (poles ⇒ ‖S₁₁‖ → ∞).
function resonance_indicator(SG::SMatrix)
    s = svdvals(SG.S11)
    log10(s[1] + 1e-15)          # σ_max
end

# ---- RCWA solver ----
# Rigorous coupled-wave analysis for a stack of laterally periodic (or homogeneous) layers
# under plane-wave illumination. The first and last entries of `layers` are the semi-infinite
# incidence and transmission half-spaces (must be HomogeneousLayer); inner layers have finite
# thickness. θ, φ are the polar/azimuthal angles of incidence (radians), and (pol_TE, pol_TM)
# are the complex TE/TM amplitudes of the incident plane wave.
function solve(layers::AbstractVector{<:AbstractLayer}, grid::RCWAGrid, f::Real;
               θ::Real=0.0, φ::Real=0.0, pol_TE=1.0, pol_TM=0.0)
    ε_ref = 𝜀̃(material(layers[1]), f)
    n_ref = sqrt(ε_ref)
    k0 = 2π * f / c₀
    kx_abs = real(n_ref) * k0 * sin(θ) * cos(φ)
    ky_abs = real(n_ref) * k0 * sin(θ) * sin(φ)
    SG, kx, ky, kz_ref, kz_trn, _, _ = smatrix_k(layers, grid, f, kx_abs, ky_abs)
    NH = nharmonics(grid)

    k̂ = (sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
    aTE = θ ≈ 0 ? (0.0, 1.0, 0.0) : (-sin(φ), cos(φ), 0.0)
    aTM = _normalize3(_cross3(aTE, k̂))
    Px = pol_TE * aTE[1] + pol_TM * aTM[1]
    Py = pol_TE * aTE[2] + pol_TM * aTM[2]

    δ = zeros(Complex{Float64}, NH)
    δ[grid.center] = 1
    csrc = vcat(Px * δ, Py * δ)

    cref = SG.S11 * csrc
    ctrn = SG.S21 * csrc
    rx = cref[1:NH]; ry = cref[(NH + 1):2NH]
    tx = ctrn[1:NH]; ty = ctrn[(NH + 1):2NH]
    kz_ref_safe = [abs(z) < 1e-12 ? Complex{Float64}(1e-12, 0) : z for z in kz_ref]
    kz_trn_safe = [abs(z) < 1e-12 ? Complex{Float64}(1e-12, 0) : z for z in kz_trn]
    rz = -Diagonal(1.0 ./ kz_ref_safe) * (Diagonal(kx) * rx + Diagonal(ky) * ry)
    tz = -Diagonal(1.0 ./ kz_trn_safe) * (Diagonal(kx) * tx + Diagonal(ky) * ty)

    kz_inc = max(real(n_ref) * cos(θ), 1e-12)
    R = [max(real(kz_ref[i]), 0.0) / kz_inc * (abs2(rx[i]) + abs2(ry[i]) + abs2(rz[i])) for i in 1:NH]
    T = [max(real(kz_trn[i]), 0.0) / kz_inc * (abs2(tx[i]) + abs2(ty[i]) + abs2(tz[i])) for i in 1:NH]
    Rtot = sum(R)
    Ttot = sum(T)
    RCWAResult(copy(grid.orders), R, T, Rtot, Ttot, 1 - Rtot - Ttot)
end


# ---- Small 3-vector helpers ----
_cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

function _normalize3(a)
    nrm = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    nrm == 0 ? a : (a[1] / nrm, a[2] / nrm, a[3] / nrm)
end
