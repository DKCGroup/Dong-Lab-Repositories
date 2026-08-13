using GESTALT
using LinearAlgebra
using JSON3
include(joinpath(@__DIR__, "metasurface_ir.jl"))
include(joinpath(@__DIR__, "config.jl"))
BLAS.set_num_threads(1)

# ---- 1. material_from_config 各类型 -----------------------------------------
println("== 1. materials ==")
for cfg in (
    Dict("type" => "nk", "n" => 1.5, "k" => 0.0),
    Dict("type" => "drude", "omega_p" => 2.4e16, "gamma_p" => 1.24e14, "eps_inf" => 1.0),
    Dict("type" => "lorentz", "eps_inf" => 8.4, "oscillators" => [Dict("deps" => 2.0, "f0" => 17.6e12, "gamma" => 1.6e12)]),
    Dict("type" => "debye", "eps_s" => 10.0, "eps_inf" => 2.0, "tau" => 1e-12),
)
    m = material_from_config(cfg)
    println("  ", cfg["type"], " => ε(10µm) = ", round(𝜀̃(m, c₀/10e-6), digits=4))
end

m = material_from_config(Dict("type" => "sampled_nk", "file" => "data/VO2-M.csv", "wavelength_unit" => "micron"))
println("  sampled_nk file => n(10µm) = ", round(real(ñ(m, c₀/10e-6)), digits=4))

# ---- 2. 完整 JSON 配置 → simulate ------------------------------------------
println("\n== 2. simulate via JSON ==")
config_json = """
{
  "period": 5e-6,
  "N": 4,
  "wavelength": {"min": 8e-6, "max": 14e-6, "n": 6},
  "incidence": {"theta": 0.0, "polarization": "unpolarized"},
  "materials": {
    "air":     {"type": "nk", "n": 1.0, "k": 0.0},
    "vo2_ins": {"type": "sampled_nk", "file": "data/VO2-I.csv", "wavelength_unit": "micron"},
    "vo2_met": {"type": "sampled_nk", "file": "data/VO2-M.csv", "wavelength_unit": "micron"},
    "al":      {"type": "drude", "omega_p": 2.4e16, "gamma_p": 1.24e14, "eps_inf": 1.0}
  },
  "layers": [
    {"material": "air", "thickness": null},
    {"material": "vo2_ins", "pattern": {"type": "shape", "shape": "circle", "thickness": 4e-6, "background": "air"}},
    {"material": "al", "thickness": 200e-9},
    {"material": "air", "thickness": null}
  ]
}
"""
out_json = simulate_json(config_json)
res = JSON3.read(out_json)
println("  polarization = ", res.polarization)
println("  n_harmonics  = ", res.n_harmonics, "   n_layers = ", res.n_layers)
println("  wavelength(µm) = ", round.(res.wavelength_um, digits=2))
println("  R = ", round.(res.R, digits=3))
println("  A = ", round.(res.A, digits=3))
println("  max_abs_RTA_deviation = ", res.max_abs_RTA_deviation)

# ---- 3. Dict 配置（多高度灰度图 + 多层介电间隔层） --------------------------
println("\n== 3. grayscale + multi-dielectric spacers ==")
img = zeros(Int, 64, 64)
for j in 1:64, i in 1:64
    r = hypot((i - 32.5) / 64, (j - 32.5) / 64)
    img[i, j] = r < 0.3 ? 200 : (r < 0.4 ? 120 : 0)   # 阶梯高度
end
config2 = Dict(
    "period" => 5e-6, "N" => 4,
    "wavelength" => Dict("min" => 8e-6, "max" => 14e-6, "n" => 5),
    "materials" => Dict(
        "air" => Dict("type" => "nk", "n" => 1.0, "k" => 0.0),
        "ge"  => Dict("type" => "nk", "n" => 4.0, "k" => 0.0),      # 锗（均匀介电间隔层）
        "vo2" => Dict("type" => "sampled_nk", "file" => "data/VO2-M.csv", "wavelength_unit" => "micron"),
    ),
    "layers" => [
        Dict("material" => "air", "thickness" => nothing),
        Dict("material" => "vo2", "pattern" => Dict("type" => "grayscale", "image" => img,
              "max_thickness" => 3e-6, "min_gray" => 50, "max_gray" => 200, "background" => "air",
              "nslices" => 4)),
        Dict("material" => "ge", "thickness" => 500e-9),             # 均匀介电间隔层
        Dict("material" => "air", "thickness" => nothing),
    ],
)
res2 = simulate(config2)
println("  n_layers = ", res2.n_layers, "  n_harmonics = ", res2.n_harmonics)
println("  A = ", round.(res2.A, digits=3))
println("  max_abs_RTA_deviation = ", res2.max_abs_RTA_deviation)

println("\nALL OK")
