# config.jl — 配置驱动的通用多层图案化超表面求解器。
#
# 将 `metasurface_ir.jl` 从「VO2 + 铝镜 + 8–14 µm」专用场景泛化为任意多层、任意材料、
# 任意波段的结构化仿真。输入是一个 JSON 配置（或等价的 Julia Dict），输出是结构化结果
# (NamedTuple)，可直接序列化为 JSON。
#
# 单位约定（全部 SI）：
#   • 长度 / 波长 / 厚度：米 (m)
#   • 频率：Hz；角频率 ωₚ、γ（Drude）：rad/s
#   • Lorentz 振子频率 f₀、线宽 γ：Hz（内部自动 ×2π）
#   • 角度：弧度；Debye 弛豫时间 τ：秒
#
# 通过 JSON 字符串跨语言调用（`simulate_json`），规避 juliacall 的类型转换细节：
#   Python  dict ──json.dumps──▶ JSON 字符串 ──▶ simulate_json ──▶ JSON 字符串 ──json.loads──▶ Python dict
#
# 见 README.md 的「通用配置驱动接口」章节了解完整 schema。

using JSON3

if !isdefined(@__MODULE__, :load_sampled_nk_csv)
    include(joinpath(@__DIR__, "metasurface_ir.jl"))
end

# Project root (for resolving relative data-file paths).
const _PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))


# ---- Generic access / conversion helpers ---------------------------------------------------

# Read a config key (String) from a Dict{Symbol,String} or a JSON3.Object; `nothing` is
# returned both for an absent key and for a JSON `null`.
function _get(cfg, key::String, default=nothing)
    cfg === nothing && return default
    haskey(cfg, key) && return cfg[key]
    haskey(cfg, Symbol(key)) && return cfg[Symbol(key)]
    default
end

# Nested JSON array (or an existing matrix) → dense Float64 matrix (for grayscale images).
function _as_matrix(x)
    x isa AbstractMatrix && return Float64.(x)
    rows = collect(x)
    nrows = length(rows)
    nrows == 0 && return zeros(Float64, 0, 0)
    row0 = collect(rows[1])
    ncols = length(row0)
    M = Matrix{Float64}(undef, nrows, ncols)
    for i in 1:nrows
        r = collect(rows[i])
        length(r) == ncols || error("image rows have inconsistent lengths (ragged array)")
        for j in 1:ncols
            M[i, j] = Float64(r[j])
        end
    end
    M
end

# Nested JSON array → Float64 vector.
_as_vector(x) = Float64.(collect(x))

# Resolve a possibly-relative data path against the project root.
_resolve_path(p::AbstractString) =
    isabspath(p) ? p : normpath(joinpath(_PROJECT_ROOT, p))


# ---- Material factory ----------------------------------------------------------------------

# Parse one material config → AbstractMaterial. Supported `type`:
#   "nk"          : {"n", "k"}
#   "drude"       : {"omega_p" [rad/s], "gamma_p" [rad/s], "eps_inf"}
#   "lorentz"     : {"eps_inf", "oscillators": [{"deps","f0"[Hz],"gamma"[Hz]}]}
#   "debye"       : {"eps_s", "eps_inf", "tau" [s]}
#   "sampled_nk"  : {"file", "wavelength_unit"} 或 内联 {"wavelength_nm","n","k"}
function material_from_config(cfg)
    typ = string(_get(cfg, "type", "nk"))
    if typ == "nk"
        return NKMaterial(Float64(_get(cfg, "n", 1.0)), Float64(_get(cfg, "k", 0.0)))
    elseif typ == "drude"
        wp = Float64(_get(cfg, "omega_p"))
        gp = Float64(_get(cfg, "gamma_p", 0.0))
        einf = Float64(_get(cfg, "eps_inf", 1.0))
        return DrudeMaterial(wp, gp, einf)
    elseif typ == "lorentz"
        einf = Float64(_get(cfg, "eps_inf"))
        oscs = _get(cfg, "oscillators", Any[])
        os = Tuple{Float64,Float64,Float64}[]
        for o in oscs
            push!(os, (Float64(_get(o, "deps")),
                       2π * Float64(_get(o, "f0")),
                       2π * Float64(_get(o, "gamma"))))
        end
        return LorentzMaterial(einf, os)
    elseif typ == "debye"
        return DebyeMaterial(Float64(_get(cfg, "eps_s")),
                             Float64(_get(cfg, "eps_inf")),
                             Float64(_get(cfg, "tau")))
    elseif typ == "sampled_nk"
        file = _get(cfg, "file", nothing)
        if file !== nothing
            unit = string(_get(cfg, "wavelength_unit", "micron"))
            λu = unit == "nm" ? 1e-9 : unit == "m" ? 1.0 : 1e-6
            return load_sampled_nk_csv(_resolve_path(String(file)); λ_unit=λu)
        else
            wl_nm = _as_vector(_get(cfg, "wavelength_nm"))
            n = _as_vector(_get(cfg, "n"))
            k = _as_vector(_get(cfg, "k"))
            return SampledDataMaterial(wl_nm, n, k; interp_kind=Gridded(Linear()))
        end
    else
        error("unknown material type '$typ'; supported: nk, drude, lorentz, debye, sampled_nk")
    end
end


# ---- Pattern factory -----------------------------------------------------------------------

# Build the patterned layer(s) from a pattern config. Returns a Vector{AbstractLayer}.
# `mat_pat` is the pattern material (from the enclosing layer's `material`); `materials` is the
# resolved name→material dict (for the `background` reference). Supported `type`:
#   "grayscale" : {"image", "max_thickness", "min_gray", "max_gray", "background", "nslices"}
#   "binary"    : {"image", "thickness", "background"}
#   "shape"     : {"shape", "thickness", "background", "pixels", "gray"}
function pattern_layers_from_config(pcfg, mat_pat::AbstractMaterial, materials::AbstractDict,
                                    period::Real, N::Integer)
    typ = string(_get(pcfg, "type", "grayscale"))
    bgname = string(_get(pcfg, "background", "air"))
    mat_bg = get(materials, bgname, get(materials, Symbol(bgname), NKMaterial(1.0)))
    mm = 2N

    if typ == "grayscale"
        img = _as_matrix(_get(pcfg, "image"))
        max_t = Float64(_get(pcfg, "max_thickness", 4e-6))
        min_gray = Float64(_get(pcfg, "min_gray", 50))
        max_gray = Float64(_get(pcfg, "max_gray", 255))
        nslices = Int(_get(pcfg, "nslices", 8))
        return build_pattern_layers(img, mat_pat, mat_bg; p=period, max_t=max_t,
                                    t_min_gray=min_gray, max_gray=max_gray, N=N,
                                    nslices=nslices)
    elseif typ == "binary"
        img = _as_matrix(_get(pcfg, "image"))
        thk = Float64(_get(pcfg, "thickness", 4e-6))
        ind = Float64.(img .> 0.0)
        pat = fourier_pattern_from_indicator(ind, period, mat_pat, mat_bg; mmax=mm, nmax=mm)
        return AbstractLayer[PeriodicLayer(mat_pat, thk, pat)]
    elseif typ == "shape"
        shp = string(_get(pcfg, "shape", "circle"))
        thk = Float64(_get(pcfg, "thickness", 4e-6))
        pixels = Int(_get(pcfg, "pixels", 64))
        gray = Int(_get(pcfg, "gray", 255))
        img = shape_image(shp; N=pixels, gray=gray)
        ind = Float64.(img .> 0.0)
        pat = fourier_pattern_from_indicator(ind, period, mat_pat, mat_bg; mmax=mm, nmax=mm)
        return AbstractLayer[PeriodicLayer(mat_pat, thk, pat)]
    else
        error("unknown pattern type '$typ'; supported: grayscale, binary, shape")
    end
end


# ---- Layer stack factory -------------------------------------------------------------------

# Parse the whole config and build (grid, layers). `layers` order: incidence half-space first,
# transmission half-space last; inner layers are homogeneous (finite thickness) or patterned.
function stack_from_config(config)
    period = Float64(_get(config, "period", 5e-6))
    N = Int(_get(config, "N", 6))

    materials_cfg = _get(config, "materials", nothing)
    materials = Dict{String,AbstractMaterial}()
    materials_cfg !== nothing && for (k, v) in pairs(materials_cfg)
        materials[string(k)] = material_from_config(v)
    end

    layers_cfg = _get(config, "layers", nothing)
    layers_cfg === nothing && error("config must contain a 'layers' array")
    nlay = length(layers_cfg)
    nlay >= 2 || error("need at least incidence and transmission half-spaces (got $nlay)")

    layers = AbstractLayer[]
    for lcfg in layers_cfg
        mname = string(_get(lcfg, "material"))
        haskey(materials, mname) || error("layer references undefined material '$mname'")
        mat = materials[mname]
        pcfg = _get(lcfg, "pattern", nothing)
        thk = _get(lcfg, "thickness", nothing)
        if pcfg !== nothing
            append!(layers, pattern_layers_from_config(pcfg, mat, materials, period, N))
        else
            d = thk === nothing ? Inf : Float64(thk)
            push!(layers, HomogeneousLayer(mat, d))
        end
    end

    grid = RCWAGrid((2π / period, 0.0), (0.0, 2π / period); Nx=N, Ny=N)
    (grid, layers)
end


# ---- Simulator -----------------------------------------------------------------------------

# Run a config dict/object and return a structured result NamedTuple. The result contains the
# wavelength grid (metres + micrometres), R/T/A for TE, TM and their combination, plus
# diagnostics (harmonics, layer count, energy-conservation deviation).
function simulate(config)
    period = Float64(_get(config, "period", 5e-6))
    grid, layers = stack_from_config(config)

    # wavelength grid: either a range {"min","max","n"} or an explicit array
    wl_cfg = _get(config, "wavelength", nothing)
    if wl_cfg !== nothing && _get(wl_cfg, "min", nothing) !== nothing
        λmin = Float64(_get(wl_cfg, "min"))
        λmax = Float64(_get(wl_cfg, "max"))
        nλ = Int(_get(wl_cfg, "n", 200))
        nλ >= 2 || error("wavelength.n must be >= 2, got $nλ")
        λs = collect(range(λmin, λmax; length=nλ))
    else
        arr = wl_cfg !== nothing ? wl_cfg : _get(config, "wavelengths", nothing)
        λs = arr === nothing ? collect(range(8e-6, 14e-6; length=200)) : _as_vector(arr)
    end

    # incidence
    inc_cfg = _get(config, "incidence", nothing)
    θ = Float64(_get(inc_cfg, "theta", 0.0))
    φ = Float64(_get(inc_cfg, "phi", 0.0))
    pol = lowercase(strip(string(_get(inc_cfg, "polarization", "unpolarized"))))

    d = sweep_spectrum_detailed(layers, grid, λs; θ=θ, φ=φ)

    if pol == "te"
        R, T, A = d.R_TE, d.T_TE, d.A_TE
    elseif pol == "tm"
        R, T, A = d.R_TM, d.T_TM, d.A_TM
    elseif pol in ("unpolarized", "unpolarised", "avg", "average")
        R = 0.5 .* (d.R_TE .+ d.R_TM)
        T = 0.5 .* (d.T_TE .+ d.T_TM)
        A = 0.5 .* (d.A_TE .+ d.A_TM)
    else
        error("unknown polarization '$pol'; use TE, TM or unpolarized")
    end

    maxdev = maximum(abs.((R .+ T .+ A) .- 1.0))

    (wavelength=λs,
     wavelength_um=λs .* 1e6,
     polarization=pol,
     R=R, T=T, A=A,
     R_TE=d.R_TE, T_TE=d.T_TE, A_TE=d.A_TE,
     R_TM=d.R_TM, T_TM=d.T_TM, A_TM=d.A_TM,
     period=period,
     n_harmonics=nharmonics(grid),
     n_layers=length(layers),
     max_abs_RTA_deviation=maxdev)
end


# ---- JSON entry points ---------------------------------------------------------------------

# Parse a JSON config string, run the simulation, and return the result as a JSON string.
function simulate_json(json_str::AbstractString)
    config = JSON3.read(json_str)
    result = simulate(config)
    JSON3.write(result)
end

# Serialize a config dict/object to a JSON string (for round-tripping / debugging).
config_to_json(config) = JSON3.write(config)

# Parse a JSON config file and return the parsed object (for Julia-side usage).
read_config_json(path::AbstractString) = JSON3.read(read(path, String))
