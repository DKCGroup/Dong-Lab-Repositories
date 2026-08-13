using LinearAlgebra

# RCWA4D — rigorous coupled-wave analysis for stacks of arbitrarily twisted periodic layers
# (any number of layers N, not only bilayers). Each device layer i carries its own reciprocal
# lattice Gᵢ and twist angle αᵢ. Following the supplement, the joint reciprocal basis is the
# direct sum G₁ ⊕ G₂ ⊕ … ⊕ G_N = { Σᵢ gᵢ : gᵢ ∈ Gᵢ }: every plane wave is labelled by one
# integer order (mᵢ, nᵢ) per layer. A layer with lattice Gᵢ only scatters within its own
# component, so its convolution matrix is block-structured (it couples two joint harmonics
# only when all other components are equal). This lets the standard RCWA machinery
# (_layer_modes, _layer_smatrix, region modes, Redheffer star product) be reused directly.
#
# The cost grows as ∏ᵢ |Gᵢ| (exponential in the number of distinct twists), as noted in the
# supplement; keep per-layer truncations small. Homogeneous spacer layers should use a
# truncation of 0 so they do not enlarge the joint basis.
if !isdefined(@__MODULE__, :AbstractLayer)
    include(joinpath(@__DIR__, "RCWA.jl"))
end


# ---- Twisted joint reciprocal grid ----
struct TwistedGrid
    """
    TwistedGrid(b1, b2, twist; Nx, Ny)
    Joint reciprocal grid for N twisted layers. All arguments are length-N vectors (one entry
    per device layer, in the same order the device layers are passed to `solve_twisted`).
    - b1, b2: Per-layer reciprocal lattice vectors (gx, gy), before twist rotation
    - twist:  Per-layer twist angle αᵢ (radians); layer i's lattice is rotated by αᵢ
    - Nx, Ny: Per-layer truncation orders (scalars are broadcast to all layers). Use 0 for
              homogeneous spacer layers so the joint basis is not enlarged.
    The joint basis enumerates every combination of per-layer orders; `gx`, `gy` hold the
    total in-plane reciprocal vector Σᵢ R(αᵢ)·(mᵢ b1ᵢ + nᵢ b2ᵢ) for each joint harmonic.
    """
    b1::Vector{NTuple{2,Float64}}
    b2::Vector{NTuple{2,Float64}}
    twist::Vector{Float64}
    multi::Vector{Vector{NTuple{2,Int}}}
    gx::Vector{Float64}
    gy::Vector{Float64}
    center::Int
    function TwistedGrid(b1, b2, twist; Nx, Ny)
        N = length(b1)
        (length(b2) == N && length(twist) == N) ||
            error("b1, b2, twist must have the same length (number of device layers)")
        Nxv = Nx isa Integer ? fill(Int(Nx), N) : collect(Int, Nx)
        Nyv = Ny isa Integer ? fill(Int(Ny), N) : collect(Int, Ny)
        (length(Nxv) == N && length(Nyv) == N) || error("Nx, Ny must be scalars or length-N")
        b1f = [NTuple{2,Float64}(b) for b in b1]
        b2f = [NTuple{2,Float64}(b) for b in b2]
        tw = collect(Float64, twist)

        per_layer = [[(m, n) for n in (-Nyv[i]):Nyv[i] for m in (-Nxv[i]):Nxv[i]] for i in 1:N]
        multi = Vector{NTuple{2,Int}}[]
        for combo in Iterators.product(per_layer...)
            push!(multi, collect(NTuple{2,Int}, combo))
        end

        gx = Float64[]
        gy = Float64[]
        center = 0
        for (p, mv) in enumerate(multi)
            tx = 0.0
            ty = 0.0
            allzero = true
            for i in 1:N
                (m, n) = mv[i]
                (m == 0 && n == 0) || (allzero = false)
                bx = m * b1f[i][1] + n * b2f[i][1]
                by = m * b1f[i][2] + n * b2f[i][2]
                c = cos(tw[i]); s = sin(tw[i])
                tx += c * bx - s * by
                ty += s * bx + c * by
            end
            push!(gx, tx)
            push!(gy, ty)
            allzero && (center = p)
        end
        center != 0 || error("zero joint order missing from grid")
        # Guard against coincident joint harmonics. At aligned (twist = 0) or commensurate twists,
        # distinct multi-indices can map to the same total g, making the joint basis G₁⊕…⊕G_N
        # degenerate (duplicated plane waves ⇒ singular convolution/eigenproblem ⇒ wrong results).
        # These cases need single-lattice RCWA (aligned) or degenerate-basis merging (not yet
        # implemented), so refuse them explicitly rather than returning garbage.
        bmax = maximum(hypot(b[1], b[2]) for b in vcat(b1f, b2f))
        tol = 1e-6 * bmax
        for p in 1:length(gx), q in (p + 1):length(gx)
            if hypot(gx[p] - gx[q], gy[p] - gy[q]) < tol
                error("TwistedGrid has coincident joint harmonics — the twist is aligned/" *
                      "commensurate and the joint basis is degenerate. Use single-lattice RCWA " *
                      "for aligned layers, or pick a generic (incommensurate) twist. " *
                      "Degenerate-basis merging for commensurate angles is not yet implemented.")
            end
        end
        new(b1f, b2f, tw, multi, gx, gy, center)
    end
end

nharmonics(grid::TwistedGrid) = length(grid.multi)
ndevicelayers(grid::TwistedGrid) = length(grid.b1)


# ---- Joint-basis convolution matrix for one device layer ----
# [εr_i] in the joint basis: layer i couples joint harmonics p and q only when their orders
# agree in every component except i; the value is then ε̂_i(Δmᵢ, Δnᵢ). Homogeneous layers
# give ε·I (diagonal in the joint basis).
function joint_convolution_matrix(layer::AbstractLayer, i::Integer, grid::TwistedGrid, f::Real)
    NH = nharmonics(grid)
    if layer isa HomogeneousLayer
        return Matrix{Complex{Float64}}(𝜀̃(material(layer), f) * I, NH, NH)
    elseif layer isa PeriodicLayer && layer.pattern isa FourierPeriodicPattern
        coeffs = _twist_order_map(layer.pattern, grid.b1[i], grid.b2[i], f)
        ER = zeros(Complex{Float64}, NH, NH)
        for p in 1:NH, q in 1:NH
            others_equal = true
            for j in 1:ndevicelayers(grid)
                j == i && continue
                grid.multi[p][j] == grid.multi[q][j] || (others_equal = false; break)
            end
            others_equal || continue
            Δm = grid.multi[p][i][1] - grid.multi[q][i][1]
            Δn = grid.multi[p][i][2] - grid.multi[q][i][2]
            ER[p, q] = get(coeffs, (Δm, Δn), 0.0 + 0.0im)
        end
        return ER
    else
        error("joint_convolution_matrix not supported for $(typeof(layer)) with this pattern")
    end
end

# Map a layer's Fourier amplitudes to integer orders (m, n) of its own (unrotated) lattice.
function _twist_order_map(pattern::FourierPeriodicPattern, b1, b2, f::Real)
    coeffs = Dict{NTuple{2,Int},Complex{Float64}}()
    B = [b1[1] b2[1]; b1[2] b2[2]]
    for (g, c) in zip(pattern.g_vectors, pattern.ε̂_coefficients)
        mn = B \ [g[1]; g[2]]
        m = round(Int, mn[1]); n = round(Int, mn[2])
        (abs(mn[1] - m) < 1e-6 && abs(mn[2] - n) < 1e-6) ||
            error("Fourier g-vector $g is not an integer combination of this layer's (b1, b2)")
        coeffs[(m, n)] = _eval_ε̂_coeff(c, f)
    end
    coeffs
end


# ---- Twisted multilayer S-matrix at fixed Bloch momentum ----
function smatrix_k_twisted(incidence::HomogeneousLayer,
                           device_layers::AbstractVector{<:AbstractLayer},
                           transmission::HomogeneousLayer,
                           grid::TwistedGrid, f::Real,
                           kx_abs::Real, ky_abs::Real)
    N = length(device_layers)
    N == ndevicelayers(grid) ||
        error("number of device layers ($N) must match grid ($(ndevicelayers(grid)))")
    NH = nharmonics(grid)
    k0 = 2π * f / c₀
    ε_ref = 𝜀̃(material(incidence), f)
    ε_trn = 𝜀̃(material(transmission), f)
    n_ref = sqrt(ε_ref)
    kx = [kx_abs / k0 + grid.gx[p] / k0 for p in 1:NH]
    ky = [ky_abs / k0 + grid.gy[p] / k0 for p in 1:NH]

    Wg, Vg, _ = _region_modes(1.0, kx, ky)
    Wref, Vref, kz_ref = _region_modes(ε_ref, kx, ky)
    Wtrn, Vtrn, kz_trn = _region_modes(ε_trn, kx, ky)

    SG = _reflection_smatrix(Wref, Vref, Wg, Vg)
    for i in 1:N
        layer = device_layers[i]
        d = thickness(layer)
        isfinite(d) || error("device layer $i must have finite thickness, got $d")
        if layer isa HomogeneousLayer
            W, V, λ = _homogeneous_modes(𝜀̃(material(layer), f), kx, ky)
        else
            ER = joint_convolution_matrix(layer, i, grid, f)
            W, V, λ = _layer_modes(ER, kx, ky)
        end
        SG = star(SG, _layer_smatrix(W, V, λ, Wg, Vg, k0, d))
    end
    SG = star(SG, _transmission_smatrix(Wtrn, Vtrn, Wg, Vg))
    (SG, kx, ky, kz_ref, kz_trn, k0, n_ref)
end


# ---- Twisted multilayer solver ----
# RCWA for a stack of arbitrarily twisted periodic (and/or homogeneous) layers. `incidence`
# and `transmission` are semi-infinite HomogeneousLayers; `device_layers` are the finite
# layers (twisted crystals and any homogeneous spacers), in order, matching `grid`. Returns an
# RCWAResult with per-(joint-)order diffraction efficiencies.
function solve_twisted(incidence::HomogeneousLayer,
                       device_layers::AbstractVector{<:AbstractLayer},
                       transmission::HomogeneousLayer,
                       grid::TwistedGrid, f::Real;
                       θ::Real=0.0, φ::Real=0.0, pol_TE=1.0, pol_TM=0.0)
    ε_ref = 𝜀̃(material(incidence), f)
    n_ref = sqrt(ε_ref)
    k0 = 2π * f / c₀
    kx_abs = real(n_ref) * k0 * sin(θ) * cos(φ)
    ky_abs = real(n_ref) * k0 * sin(θ) * sin(φ)
    SG, kx, ky, kz_ref, kz_trn, _, _ =
        smatrix_k_twisted(incidence, device_layers, transmission, grid, f, kx_abs, ky_abs)
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
    R = [max(real(kz_ref[p]), 0.0) / kz_inc * (abs2(rx[p]) + abs2(ry[p]) + abs2(rz[p])) for p in 1:NH]
    T = [max(real(kz_trn[p]), 0.0) / kz_inc * (abs2(tx[p]) + abs2(ty[p]) + abs2(tz[p])) for p in 1:NH]
    Rtot = sum(R); Ttot = sum(T)
    orders = [grid.multi[p][1] for p in 1:NH]
    RCWAResult(orders, R, T, Rtot, Ttot, 1 - Rtot - Ttot)
end
