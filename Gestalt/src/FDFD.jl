# FDFD.jl — 实空间（像素网格）面内本征模求解器。
#
# 动机：RCWA 把整层 ε(r) 做全局 Fourier 展开，表示宽度 w 的特征需要谐波阶 N ≈ p/w，
# 对"大周期 + 细结构"（细介质条、细金属条）代价按 O(NH³) 爆炸且 Gibbs 截断误差不收敛。
# 本文件提供一条与 Fourier RCWA 并列的路径：在像素网格上离散 ε（精确对角）与导数
# （稀疏有限差分），用稀疏迭代本征（Arpack 位移求逆）求"最不易衰逝"的面内本征模，
# 再通过最小二乘模匹配接入现有的 S 矩阵级联。复杂度从稠密 O(NH³) 降为稀疏 ~线性。
#
# 关键接口（与现有 `_layer_modes` / `_layer_smatrix` 完全对齐）：
#   layer_modes_realspace(...) -> (W, V, λ)     实空间基（2Npix² × 2NH）
#   gap_modes_realspace(...)   -> (Wg, Vg)      间隙平面波采样到像素网格（2Npix² × 2NH）
# 二者同维后，直接复用 RCWA.jl 的 `_layer_smatrix`（其对高瘦矩阵做最小二乘 W\Wg）。

using SparseArrays
using Arpack
using LinearAlgebra
using FFTW

# ---- 1D 周期性有限差分导数矩阵 --------------------------------------------------
# 中心差分，order ∈ {2,4,6,8}。对周期格点 f_i（i=1..Npix）返回 ∂/∂x 的稀疏矩阵，
# 满足 Dx·e^{i g x} ≈ i g·e^{i g x}（即 Dx 近似微分算子 ∂/∂x）。
function fd_derivative_1d(Npix::Integer, dx::Real; order::Integer=4)
    c = order == 2 ? [0.5] :
        order == 4 ? [2 / 3, -1 / 12] :
        order == 6 ? [3 / 4, -3 / 20, 1 / 60] :
        order == 8 ? [4 / 5, -1 / 5, 4 / 105, -1 / 280] :
        error("order 必须是 2/4/6/8")
    m = length(c)
    rows = Int[]; cols = Int[]; vals = Float64[]
    for i in 1:Npix, k in 1:m
        push!(rows, i); push!(cols, mod1(i + k, Npix)); push!(vals, c[k] / dx)
        push!(rows, i); push!(cols, mod1(i - k, Npix)); push!(vals, -c[k] / dx)
    end
    sparse(rows, cols, vals, Npix, Npix)
end

# ---- 1D 周期性谱（DFT）导数矩阵 --------------------------------------------------
# 谱导数对所有空间频率都精确（无 FD 色散误差），代价是稠密矩阵。
# 对周期格点 f_i（i=1..Npix），D = F⁻¹·diag(i·g)·F，满足 D·e^{i g x} = i g·e^{i g x}。
function spectral_derivative_1d(Npix::Integer, dx::Real)
    N = Int(Npix)
    # DFT 频率（匹配 FFT 的 0..N-1 折叠顺序）
    g = [2π * (m <= N ÷ 2 ? m : m - N) / (N * dx) for m in 0:(N - 1)]
    F = [exp(-2π * im * (m - 1) * (n - 1) / N) for m in 1:N, n in 1:N]  # DFT 矩阵
    D = (1.0 / N) * (F' * Diagonal(im .* g) * F)                        # F⁻¹ = F'/N
    Matrix{Complex{Float64}}(D)
end

# ---- 谱导数 + 无矩阵 P/Q 应用（用于 Rayleigh 商修正，可扩展到 Npix=256+） ---------
# 谱导数对所有频率精确但稠密（O(Npix⁴) 存储）；这里用 FFT 无矩阵施加 ∂/∂x、∂/∂y，
# 使"稀疏 FD 本征 + 谱修正本征值"的路径能扩展到大规模像素网格。
_freq_vector(N::Integer, dx::Real) =
    [2π * (m <= N ÷ 2 ? m : m - N) / (N * dx) for m in 0:(N - 1)]

# 谱 ∂/∂x、∂/∂y：对 Npix×Npix 的 2D 场（x=第 1 维快索引，y=第 2 维）做 FFT 微分。
_spec_Dx(V::AbstractMatrix, gx::AbstractVector) =
    ifft(fft(V, 1) .* reshape(im .* gx, :, 1), 1)
_spec_Dy(V::AbstractMatrix, gy::AbstractVector) =
    ifft(fft(V, 2) .* reshape(im .* gy, 1, :), 2)

# 对扁平 2Npix² 矢量 v = [vx; vy] 施加 Q 算子（谱导数，无矩阵）。
# 与 realspace_PQ 中 Q = [KX·KY, ER−KX·KX; KY·KY−ER, −KY·KX] 一致。
function _apply_Q_spec(v::AbstractVector{<:Complex}, εvec, k0, kx_abs, ky_abs, p, Npix)
    N2 = Int(Npix^2)
    g = _freq_vector(Int(Npix), p / Npix)
    vx = copy(reshape(v[1:N2], Npix, Npix))
    vy = copy(reshape(v[N2+1:2N2], Npix, Npix))
    cx = kx_abs / k0; cy = ky_abs / k0; c = -im / k0
    ε = reshape(εvec, Npix, Npix)
    KXvx = cx .* vx .+ c .* _spec_Dx(vx, g)
    KYvx = cy .* vx .+ c .* _spec_Dy(vx, g)
    KXvy = cx .* vy .+ c .* _spec_Dx(vy, g)
    KYvy = cy .* vy .+ c .* _spec_Dy(vy, g)
    KXKYvx = cx .* KYvx .+ c .* _spec_Dx(KYvx, g)
    KXKXvy = cx .* KXvy .+ c .* _spec_Dx(KXvy, g)
    qx = KXKYvx .+ (ε .* vy .- KXKXvy)
    KYKYvx = cy .* KYvx .+ c .* _spec_Dy(KYvx, g)
    KYKXvy = cy .* KXvy .+ c .* _spec_Dy(KXvy, g)
    qy = (KYKYvx .- ε .* vx) .- KYKXvy
    vcat(vec(qx), vec(qy))
end

# 对扁平 2Npix² 矢量施加 P 算子（谱导数，无矩阵）。
# P = [KX·ERi·KY, I − KX·ERi·KX; KY·ERi·KY − I, −KY·ERi·KX]
function _apply_P_spec(v::AbstractVector{<:Complex}, εvec, k0, kx_abs, ky_abs, p, Npix)
    N2 = Int(Npix^2)
    g = _freq_vector(Int(Npix), p / Npix)
    vx = copy(reshape(v[1:N2], Npix, Npix))
    vy = copy(reshape(v[N2+1:2N2], Npix, Npix))
    cx = kx_abs / k0; cy = ky_abs / k0; c = -im / k0
    ε = reshape(εvec, Npix, Npix)
    εi = 1.0 ./ ε
    KXvx = cx .* vx .+ c .* _spec_Dx(vx, g)
    KYvx = cy .* vx .+ c .* _spec_Dy(vx, g)
    KXvy = cx .* vy .+ c .* _spec_Dx(vy, g)
    KYvy = cy .* vy .+ c .* _spec_Dy(vy, g)
    ERiKYvx = εi .* KYvx
    KX_ERiKYvx = cx .* ERiKYvx .+ c .* _spec_Dx(ERiKYvx, g)
    ERiKXvy = εi .* KXvy
    KX_ERiKXvy = cx .* ERiKXvy .+ c .* _spec_Dx(ERiKXvy, g)
    px = KX_ERiKYvx .+ (vy .- KX_ERiKXvy)
    KY_ERiKYvx = cy .* ERiKYvx .+ c .* _spec_Dy(ERiKYvx, g)
    KY_ERiKXvy = cy .* ERiKXvy .+ c .* _spec_Dy(ERiKXvy, g)
    py = (KY_ERiKYvx .- vx) .- KY_ERiKXvy
    vcat(vec(px), vec(py))
end

# ---- 实空间 P, Q 算子 -----------------------------------------------------------
# 与 RCWA.jl `_layer_modes` 中稠密 P, Q 完全相同的分块结构，只是把：
#   KX, KY（归一化面内导数算子）由对角 Fourier 波矢换成 (kx_abs/k0)·I + (-i/k0)·Dx/y；
#   ER（ε 卷积矩阵）换成像素点乘的对角阵；ERi 换成逐点倒数对角阵。
# `deriv ∈ {:spectral, :fd}`：谱导数精确但稠密；FD 稀疏但高频有色散误差（用于大规模）。
function realspace_PQ(εvec::AbstractVector{<:Complex}, k0::Real,
                      kx_abs::Real, ky_abs::Real, p::Real, Npix::Integer;
                      deriv::Symbol=:spectral, order::Integer=8)
    dx = p / Npix
    if deriv == :spectral
        D1 = spectral_derivative_1d(Npix, dx)
        I1 = Matrix{Complex{Float64}}(I, Npix, Npix)
        Dx = kron(I1, D1)      # ∂/∂x
        Dy = kron(D1, I1)      # ∂/∂y
        Ibig = Matrix{Complex{Float64}}(I, Npix^2, Npix^2)
    else
        D1 = fd_derivative_1d(Npix, dx; order=order)
        I1 = sparse(1.0 * I, Npix, Npix)
        Dx = kron(I1, D1)
        Dy = kron(D1, I1)
        Ibig = sparse(1.0 * I, Npix^2, Npix^2)
    end
    KX = (-1im / k0) * Dx + (kx_abs / k0) * Ibig
    KY = (-1im / k0) * Dy + (ky_abs / k0) * Ibig
    if deriv == :spectral
        ER = Diagonal(Complex{Float64}.(εvec))
        ERi = Diagonal(1.0 ./ Complex{Float64}.(εvec))
    else
        ER = spdiagm(0 => Complex{Float64}.(εvec))
        ERi = spdiagm(0 => 1.0 ./ Complex{Float64}.(εvec))
    end

    P = [KX * ERi * KY        Ibig - KX * ERi * KX;
         KY * ERi * KY - Ibig -KY * ERi * KX]
    Q = [KX * KY              ER - KX * KX;
         KY * KY - ER         -KY * KX]
    (P, Q)
end

# ---- 实空间面内本征模 ------------------------------------------------------------
# 求解 P·Q 的最不易衰逝（最小 |λ|）本征对，返回 (W, V, λ) 全在实空间像素基：
#   W: 2Npix² × nmodes 电模；V = Q·W·diag(1/λ) 磁模；λ 取与稠密代码一致的"全前向"分支。
function layer_modes_realspace(εvec::AbstractVector{<:Complex}, k0::Real,
                               kx_abs::Real, ky_abs::Real, p::Real, Npix::Integer;
                               nmodes::Integer, deriv::Symbol=:spectral, order::Integer=8,
                               sigma::Real=1e-6, tol::Real=1e-8, maxiter::Integer=2000,
                               over::Integer=3)
    P, Q = realspace_PQ(εvec, k0, kx_abs, ky_abs, p, Npix; deriv=deriv, order=order)
    PQ = P * Q
    # 过采样：Arpack 位移求逆找的是最小 |λ²| 模，会把衰逝模（λ²>0 小值）排在传播模
    # （λ²<0）之前；为保证所有传播模（Re λ≈0）都被覆盖，多取 over·nmodes 个再按
    # 衰减率 Re λ 重排并截断。
    n = min(over * nmodes, size(PQ, 1) - 1)
    λsq, W = eigs(PQ; nev=n, which=:LM, sigma=sigma, tol=tol, maxiter=maxiter)
    λ = sqrt.(Complex{Float64}.(λsq))
    _branch_tol = 1e-8
    _branch!(λ) = for i in eachindex(λ)
        r = real(λ[i]); im = imag(λ[i])
        (r < -_branch_tol || (abs(r) <= _branch_tol && im > 0)) && (λ[i] = -λ[i])
    end
    _branch!(λ)
    if deriv == :fd
        # Rayleigh 商修正：FD 导数在高频有色散误差（含 Nyquist 处 D̂=0 的虚假"传播"模）。
        # 用谱算子（FFT，精确）重算每个模的 λ²，把虚假模修正回其真实的大衰减率，
        # 从而在下面的 Re λ 排序中被正确放到衰逝区（不被选中）。
        for i in eachindex(λ)
            w = W[:, i]
            a = _apply_P_spec(_apply_Q_spec(w, εvec, k0, kx_abs, ky_abs, p, Npix),
                              εvec, k0, kx_abs, ky_abs, p, Npix)
            λ[i] = sqrt(Complex(dot(w, a) / dot(w, w)))
        end
        _branch!(λ)
    end
    # 按衰减率 Re λ 升序（最不易衰逝在前），传播模（Re λ≈0）优先于衰逝模；
    # 并列（同为传播模）时按 |λ| 升序。不能用 |λ| 排序——衰逝模（λ 纯实）可能比
    # 传播模（λ 纯虚）有更小的模，导致传播模被挤掉、基矢量缺失（均匀层 N=1 即源于此）。
    perm = sortperm(eachindex(λ); by=i -> (real(λ[i]), abs(λ[i])))
    λ = λ[perm]
    W = W[:, perm]
    nkeep = min(nmodes, length(λ))
    λ = λ[1:nkeep]
    W = W[:, 1:nkeep]
    if deriv == :fd
        # V 用谱 Q 计算，与修正后的 λ 一致（FD 基的 Q 会残留色散误差）。
        V = zeros(Complex{Float64}, 2 * Npix^2, nkeep)
        for i in 1:nkeep
            V[:, i] = _apply_Q_spec(W[:, i], εvec, k0, kx_abs, ky_abs, p, Npix) / λ[i]
        end
    else
        V = Q * W * Diagonal(1.0 ./ λ)
    end
    (W, V, λ)
end

# ---- 间隙平面波模采样到像素网格 ------------------------------------------------
# 间隙（ε=1 半空间）本征模在平面波基下为 Wg=I、Vg=Q·λ⁻¹。这里把平面波
# （周期部分 e^{i g·r}）采样到像素网格，得到 2Npix² × 2NH 的实空间基表示。
# 返回 (Wg_rs, Vg_rs)。
function gap_modes_realspace(grid, k0::Real, kx_abs::Real, ky_abs::Real,
                             p::Real, Npix::Integer)
    NH = nharmonics(grid)
    kx, ky = _kxky_from_kabs(grid, k0, kx_abs, ky_abs)
    _, Vg_pw, _ = _region_modes(1.0, kx, ky)          # Vg 平面波基 (2NH × 2NH)
    dx = p / Npix
    S = zeros(Complex{Float64}, 2 * Npix^2, 2 * NH)    # 采样矩阵（电模 Wg=I 的采样）
    N2 = Npix^2
    @inbounds for j in 1:Npix
        y = (j - 0.5) * dx
        base = (j - 1) * Npix
        for i in 1:Npix
            x = (i - 0.5) * dx
            q = base + i
            for h in 1:NH
                ph = cis(grid.gx[h] * x + grid.gy[h] * y)
                S[q, h] = ph                 # Ex（x 偏振谐波）
                S[q + N2, NH + h] = ph       # Ey（y 偏振谐波）
            end
        end
    end
    (S, S * Vg_pw)
end

# ---- 实空间图案层的完整 S 矩阵（间隙基，2NH × 2NH） ------------------------------
# indicator: Npix × Npix 指示矩阵（∈[0,1]，1=材料 in）；eps_in/eps_out 为当前频率下
# 两种材料的复介电常数。返回可直接并入 `smatrix_k` 级联的 SMatrix。
function layer_smatrix_realspace(indicator::AbstractMatrix, eps_in::Complex, eps_out::Complex,
                                 grid, k0::Real, kx_abs::Real, ky_abs::Real,
                                 d::Real, p::Real, Npix::Integer;
                                 deriv::Symbol=:spectral, order::Integer=8, sigma::Real=1e-6)
    NH = nharmonics(grid)
    εmat = eps_out .+ (eps_in - eps_out) .* Float64.(indicator)
    εvec = vec(Complex{Float64}.(εmat))
    W, V, λ = layer_modes_realspace(εvec, k0, kx_abs, ky_abs, p, Npix;
                                    nmodes=2 * NH, deriv=deriv, order=order, sigma=sigma)
    Wg, Vg = gap_modes_realspace(grid, k0, kx_abs, ky_abs, p, Npix)
    _layer_smatrix(W, V, λ, Wg, Vg, k0, d)
end
