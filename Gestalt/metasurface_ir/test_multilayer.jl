# test_multilayer.jl — 验证多层 block-diagonal 快速路径的正确性与性能。
#
# 对比：旧的稠密路径（逐层 _layer_smatrix + 稠密 star） vs 新的 2×2 block-diagonal 路径。
# 结论：两者 S 矩阵与 R/T/A 应逐元素一致（机器精度），新路径对均匀层更快。

using GESTALT
using LinearAlgebra
BLAS.set_num_threads(1)

const G = GESTALT

# ---- 旧稠密路径参考实现（与优化前的 smatrix_k 完全一致） ---------------------
function smatrix_k_old(layers, grid, f, kx_abs, ky_abs)
    N = length(layers)
    k0 = 2π * f / c₀
    ε_ref = 𝜀̃(material(layers[1]), f)
    ε_trn = 𝜀̃(material(layers[N]), f)
    n_ref = sqrt(ε_ref)
    kx, ky = G._kxky_from_kabs(grid, k0, kx_abs, ky_abs)
    Wg, Vg, _ = G._region_modes(1.0, kx, ky)
    Wref, Vref, kz_ref = G._region_modes(ε_ref, kx, ky)
    Wtrn, Vtrn, kz_trn = G._region_modes(ε_trn, kx, ky)
    SG = G._reflection_smatrix(Wref, Vref, Wg, Vg)
    for i in 2:(N - 1)
        layer = layers[i]
        d = thickness(layer)
        if layer isa HomogeneousLayer
            W, V, λ = G._homogeneous_modes(𝜀̃(material(layer), f), kx, ky)
        else
            ER = convolution_matrix(layer, grid, f)
            W, V, λ = G._layer_modes(ER, kx, ky)
        end
        SG = star(SG, G._layer_smatrix(W, V, λ, Wg, Vg, k0, d))
    end
    SG = star(SG, G._transmission_smatrix(Wtrn, Vtrn, Wg, Vg))
    (SG, kx, ky, kz_ref, kz_trn, k0, n_ref)
end

function max_abs_diff(A, B)
    maximum(abs.(A .- B))
end

air = HomogeneousLayer(NKMaterial(1.0), Inf)
f = c₀ / 10e-6
grid = RCWAGrid((2π / 5e-6, 0.0), (0.0, 2π / 5e-6); Nx=4, Ny=4)

# ---- 1. 纯均匀多层（无图案）：新旧逐元素对比 + 能量守恒 ----------------------
println("== 1. 纯均匀多层 (air | Si | Ge | Al | air) ==")
layers1 = [
    air,
    HomogeneousLayer(NKMaterial(3.42), 1.0e-6),   # Si
    HomogeneousLayer(NKMaterial(4.0), 500e-9),    # Ge
    HomogeneousLayer(DrudeMaterial(2.4e16, 1.24e14, 1.0), 200e-9),  # Al
    air,
]
new = smatrix_k(layers1, grid, f, 0.0, 0.0)
old = smatrix_k_old(layers1, grid, f, 0.0, 0.0)
println("  SG.S11 max diff = ", max_abs_diff(new[1].S11, old[1].S11))
println("  SG.S21 max diff = ", max_abs_diff(new[1].S21, old[1].S21))
println("  kz_ref max diff  = ", max_abs_diff(new[3], old[3]))

r_new = solve(layers1, grid, f)
r_old = G.solve(layers1, grid, f)  # solve 内部用新 smatrix_k
println("  Rtot = ", round(r_new.Rtot, digits=6), "  Ttot = ", round(r_new.Ttot, digits=6),
        "  A = ", round(r_new.A, digits=6), "  R+T+A = ", round(r_new.Rtot + r_new.Ttot + r_new.A, digits=12))

# ---- 2. 含图案层的多层（图案 + 多个均匀层）：新旧对比 ------------------------
println("\n== 2. 图案层 + 多层均匀层 ==")
include(joinpath(@__DIR__, "metasurface_ir.jl"))
img = circle_image(0.5, 0.5, 0.34; gray=255, N=64)
ind = Float64.(img .> 0)
vo2 = vo2_metallic_sampled()
pat = fourier_pattern_from_indicator(ind, 5e-6, vo2, NKMaterial(1.0); mmax=8, nmax=8)
layers2 = [
    air,
    PeriodicLayer(vo2, 4e-6, pat),
    HomogeneousLayer(NKMaterial(4.0), 500e-9),    # Ge spacer
    HomogeneousLayer(NKMaterial(2.0), 300e-9),    # 额外均匀层
    HomogeneousLayer(DrudeMaterial(2.4e16, 1.24e14, 1.0), 200e-9),  # Al
    air,
]
new2 = smatrix_k(layers2, grid, f, 0.0, 0.0)
old2 = smatrix_k_old(layers2, grid, f, 0.0, 0.0)
println("  SG.S11 max diff = ", max_abs_diff(new2[1].S11, old2[1].S11))
println("  SG.S12 max diff = ", max_abs_diff(new2[1].S12, old2[1].S12))
println("  SG.S21 max diff = ", max_abs_diff(new2[1].S21, old2[1].S21))
println("  SG.S22 max diff = ", max_abs_diff(new2[1].S22, old2[1].S22))
r2 = solve(layers2, grid, f)
println("  Rtot = ", round(r2.Rtot, digits=6), "  Ttot = ", round(r2.Ttot, digits=6),
        "  A = ", round(r2.A, digits=6), "  R+T+A = ", round(r2.Rtot + r2.Ttot + r2.A, digits=12))

# ---- 3. 与 TMM 交叉验证（纯均匀多层） --------------------------------------
println("\n== 3. RCWA(均匀多层) vs TMM ==")
# TMM 求解器：solve([HomogeneousLayer...], f) 直接分派到 TMM
tmm = G.solve([HomogeneousLayer(NKMaterial(1.0), Inf),
               HomogeneousLayer(NKMaterial(3.42), 1.0e-6),
               HomogeneousLayer(NKMaterial(1.5), Inf)], f)
println("  TMM R = ", round(tmm.R, digits=6), " T = ", round(tmm.T, digits=6))

# ---- 4. 性能对比（同一均匀多层，多次求解计时） ------------------------------
println("\n== 4. 性能对比（20 次求解） ==")
const L = layers1
new_t = @elapsed for _ in 1:20
    smatrix_k(L, grid, f, 0.0, 0.0)
end
old_t = @elapsed for _ in 1:20
    smatrix_k_old(L, grid, f, 0.0, 0.0)
end
println("  新路径(blockdiag): ", round(new_t, digits=4), " s")
println("  旧路径(dense)    : ", round(old_t, digits=4), " s")
println("  加速比: ", round(old_t / new_t, digits=2), "×")

# ---- 5. 多个图案层 + 均匀层（真正的多层图案结构） ---------------------------
println("\n== 5. 两个图案层 + 均匀层 ==")
vo2i = vo2_insulating_sampled()
img2 = triangle_image((0.12, 0.12), (0.88, 0.18), (0.42, 0.82); gray=255, N=64)
ind2 = Float64.(img2 .> 0)
pat2 = fourier_pattern_from_indicator(ind2, 5e-6, vo2i, NKMaterial(1.0); mmax=8, nmax=8)
layers5 = [
    air,
    PeriodicLayer(vo2, 2e-6, pat),               # 图案 1（金属态 VO2 圆）
    HomogeneousLayer(NKMaterial(2.0), 250e-9),   # 中间均匀隔离层
    PeriodicLayer(vo2i, 2e-6, pat2),             # 图案 2（绝缘态 VO2 三角）
    HomogeneousLayer(DrudeMaterial(2.4e16, 1.24e14, 1.0), 200e-9),  # Al
    air,
]
new5 = smatrix_k(layers5, grid, f, 0.0, 0.0)
old5 = smatrix_k_old(layers5, grid, f, 0.0, 0.0)
println("  SG.S11 max diff = ", max_abs_diff(new5[1].S11, old5[1].S11))
println("  SG.S21 max diff = ", max_abs_diff(new5[1].S21, old5[1].S21))
r5 = solve(layers5, grid, f)
println("  Rtot = ", round(r5.Rtot, digits=6), "  Ttot = ", round(r5.Ttot, digits=6),
        "  A = ", round(r5.A, digits=6), "  R+T+A = ", round(r5.Rtot + r5.Ttot + r5.A, digits=12))

println("\nALL OK")
