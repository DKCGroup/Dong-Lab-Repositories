# metasurface_ir 归档说明（维护者文档）

> 本文件是 **metasurface_ir** 子模块的归档说明，供后续回到本模块修改时快速恢复上下文，
> 无需通读全部代码。面向**使用者**的接口/JSON schema/性能说明见同目录
> [`README.md`](README.md)；本文件聚焦**代码结构、数据流、依赖关系、历史坑与验证状态**。

---

## 1. 这是什么

`metasurface_ir/` 是 GESTALT 框架下、面向 **8–14 µm 长波红外（LWIR）超表面** 的仿真工作流，
构建在核心 RCWA 引擎（`src/RCWA.jl`）之上。核心能力：

- 把**灰度图**（N×N，0–255）映射为顶层图案的**厚度分布**（0–255 → 0–`max_t`，≤阈值视为空气）。
- 二值图 → 单层 `PeriodicLayer`（快路径）；多高度灰度图 → `nslices` 个等厚阶梯子层（staircase）。
- 标准两层结构 `air | VO2 图案 | (可选均匀介电间隔层) | 铝镜(200nm) | air`；铝镜在 LWIR 不透明，故 `A = 1 − R`。
- 集成实测 VO₂ n,k 数据（`data/VO2-I.csv` 绝缘态、`data/VO2-M.csv` 金属态）。
- 提供 **Julia** 与 **Python（juliacall）** 两套调用接口；并通过 **JSON 配置** 泛化为任意多层/任意材料/任意波段。

---

## 2. 目录结构与文件职责

```
GESTALT/
├── src/                          # 核心引擎（GESTALT 模块，metasurface_ir 依赖它）
│   ├── GESTALT.jl                # 模块入口：include 顺序 + 全部 export
│   ├── Material.jl               # 材料模型（AbstractMaterial 及 7 个子类型）
│   ├── Layer.jl                  # 层类型（Homogeneous / Periodic / Aperiodic + 图案）
│   ├── TMM.jl                    # 均匀多层传输矩阵法（solve 的同名重载）
│   ├── RCWA.jl                   # 严格耦合波分析（RCWAGrid / SMatrix / star / solve）
│   ├── RCWA4D.jl                 # N 层扭转 RCWA（本模块暂未使用）
│   ├── CMT.jl / Visualize.jl     # 耦合模理论 / 可视化（本模块未直接使用）
├── data/
│   ├── VO2-I.csv                 # 绝缘态 VO2 n,k（列：波长µm, n, k；降序）
│   └── VO2-M.csv                 # 金属态 VO2 n,k（同上）
└── metasurface_ir/
    ├── metasurface_ir.jl         # ★ 核心库（材料工厂/灰度→图案/堆栈/扫频/形状）
    ├── api.jl                    # 跨语言便捷 API（metasurface_spectrum 等）
    ├── config.jl                 # ★ 配置驱动通用求解器（JSON schema → simulate）
    ├── validate_rcwa.jl          # RCWA 正确性验证（能量守恒/TMM 交叉/收敛）
    ├── test_config.jl            # 配置接口回归测试（材料/图案/JSON 往返）
    ├── test_multilayer.jl        # 多层 block-diagonal 快速路径正确性 + 性能
    ├── benchmark.jl              # 批量生成性能基准（N/n_lambda/nslices 耗时）
    ├── run_shapes.jl             # 形状仿真 + 可视化（矩形/圆/三角/十字）
    ├── figures/                  # 输出图像（含 beamer 用图）
    ├── README.md                 # 使用者文档（JSON schema / API / 性能）
    ├── ARCHIVE.md                # 本文件
    └── python/
        ├── gestalt_metasurface.py  # ★ Python API 封装（juliacall）
        ├── example.py              # Python 使用示例
        ├── test_config.py          # Python 端配置接口回归测试
        ├── make_beamer_figures.py  # 生成 beamer 演示图
        ├── demo.ipynb              # Jupyter 演示 notebook
        └── requirements.txt        # juliacall / numpy / matplotlib / jupyter
```

> 标记 ★ 的是理解本模块时**必读**的文件（共 3 个：`metasurface_ir.jl`、`config.jl`、`gestalt_metasurface.py`）。

---

## 3. 依赖关系与 include 顺序

三个 Julia 文件都是**脚本（include 使用）而非子模块**，函数进入调用方的全局作用域：

```
metasurface_ir.jl   依赖  src 全部（需先 `using GESTALT`）
api.jl              依赖  metasurface_ir.jl（文件内自含 include 守卫）
config.jl           依赖  metasurface_ir.jl（文件内自含 include 守卫）+ JSON3
```

- 每个文件顶部都有 `if !isdefined(@__MODULE__, :xxx); include(...); end` 守卫，可安全重复 include。
- `metasurface_ir.jl` 内 `const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))`，
  `config.jl` 内 `const _PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))` 用于解析相对路径。
- 核心引擎的类型/函数在 `src/GESTALT.jl` 全部 export，`metasurface_ir` 直接使用，**不要重定义**。

---

## 4. 两条调用路径（数据流）

### 路径 A：便捷 API（VO2 + 铝镜专用）

```
Python: run_shape("circle", "metallic")  /  run_spectrum(img, phase)
   │  juliacall（关键字参数名必须 ASCII）
   ▼
Julia api.jl: metasurface_shape_spectrum(shape, phase; kwargs)
   └─ metasurface_spectrum(shape_image(shape), phase; kwargs)
        └─ _phase_material(phase; sampled)  →  vo2_{insulating|metallic}_sampled/material()
        └─ build_stack(img, mat, aluminum_material(); p, max_t, N, d_al)
             └─ build_pattern_layers(img, mat, NKMaterial(1.0); ...)
                  └─ fourier_pattern_from_indicator(ind, p, mat_in, mat_out; mmax=2N, nmax=2N)
        └─ sweep_spectrum(layers, grid, λs; θ)   # Threads.@threads 并行
             └─ solve(layers, grid, f; θ, pol_TE, pol_TM)   # src/RCWA.jl
```

### 路径 B：通用 JSON 配置接口

```
Python: simulate(config_dict)
   │  json.dumps(config)  ──▶  JSON 字符串
   ▼
Julia config.jl: simulate_json(json_str)
   └─ simulate(JSON3.read(json_str))
        ├─ stack_from_config(config)
        │    ├─ material_from_config(cfg)         # nk/drude/lorentz/debye/sampled_nk
        │    ├─ pattern_layers_from_config(pcfg, ...)  # grayscale/binary/shape
        │    └─ RCWAGrid(...; Nx=N, Ny=N)
        ├─ sweep_spectrum_detailed(layers, grid, λs; θ, φ)   # 返回 TE/TM 分离
        └─ 按 polarization 组合 → NamedTuple 结果
   │  JSON3.write(result)  ──▶  JSON 字符串
   ▼
Python: json.loads(...)  ──▶  dict
```

> 跨语言只传 **JSON 字符串**，规避 juliacall 的 NamedTuple/Dict 类型转换细节。这是设计要点，
> 修改时请勿绕过 `simulate_json` 直接返回 Julia NamedTuple 给 Python。

---

## 5. 关键 API 清单

### 5.1 metasurface_ir.jl（核心库）

| 函数 | 作用 |
|---|---|
| `vo2_insulating_material()` | Lorentz 单振子绝缘态 VO2（`LorentzMaterial(8.4, [(2.0, 2π·17.6e12, 2π·1.6e12)])`） |
| `vo2_metallic_material()` | Drude 金属态 VO2（`DrudeMaterial(1.51e15, 5.75e14, 9.0)`） |
| `aluminum_material()` | 铝镜 Drude（`DrudeMaterial(2.4e16, 1.24e14, 1.0)`） |
| `vo2_insulating_sampled()` / `vo2_metallic_sampled()` | 实测 CSV n,k（`SampledDataMaterial`） |
| `load_sampled_nk_csv(path; λ_unit=1e-6, interp_kind)` | CSV→材料；自动降序→升序、µm→nm、去重 |
| `gray_to_thickness(gray, max_t, t_min_gray, max_gray)` | 单灰度值→厚度 |
| `thickness_profile(img, max_t, t_min_gray, max_gray)` | 灰度图→厚度矩阵 |
| `fourier_pattern_from_indicator(ind, p, mat_in, mat_out; mmax, nmax)` | 二值指示函数→`FourierPeriodicPattern` |
| `build_pattern_layers(img, mat_vo2, mat_bg; p, max_t, t_min_gray, max_gray, N, nslices)` | 图像→图案层（二值单层 / 灰度阶梯） |
| `build_stack(img, mat_vo2, mat_al; ...)` | 图像→完整层堆栈 `(grid, layers)` |
| `sweep_spectrum(layers, grid, λs; θ, φ, unpolarized)` | R/T/A 扫频（unpolarized 时 TE/TM 平均） |
| `sweep_spectrum_detailed(layers, grid, λs; θ, φ)` | TE/TM 分离扫频 |
| `shape_image(shape; N, gray)` / `rect_image` / `circle_image` / `triangle_image` / `cross_image` | 内置形状灰度图生成器 |

### 5.2 api.jl（便捷封装）

| 函数 | 作用 |
|---|---|
| `metasurface_spectrum(img, phase; p, max_t, t_min_gray, N, lambda_min, lambda_max, n_lambda, theta, sampled, d_al)` | 灰度图→R/T/A |
| `metasurface_shape_spectrum(shape, phase; kwargs...)` | 形状名→R/T/A |
| `shape_image(shape; N, gray)` | 形状→灰度图（转发到核心库） |
| `vo2_dispersion(phase; lambdas, sampled)` | 复折射率 n+ik |

### 5.3 config.jl（通用接口）

| 函数 | 作用 |
|---|---|
| `material_from_config(cfg)` | 材料定义 dict → `AbstractMaterial` |
| `pattern_layers_from_config(pcfg, mat_pat, materials, period, N)` | 图案定义 → `Vector{AbstractLayer}` |
| `stack_from_config(config)` | 配置 → `(grid, layers)` |
| `simulate(config)` | 配置 → 结构化结果 `NamedTuple` |
| `simulate_json(json_str)` | JSON 字符串 → JSON 字符串（跨语言主入口） |
| `config_to_json(config)` / `read_config_json(path)` | 序列化 / 读文件 |

### 5.4 python/gestalt_metasurface.py

| 函数 | 作用 |
|---|---|
| `run_spectrum(img, phase, *, p, max_t, ..., N, n_lambda, ...)` | 灰度图→dict |
| `run_shape(shape, phase, **kwargs)` | 形状→dict |
| `shape_image(shape, N=64, gray=255)` | 形状→numpy 数组 |
| `vo2_dispersion(phase, ...)` | 色散→dict |
| `simulate(config)` / `simulate_from_json(path)` / `save_result(result, path)` | 通用 JSON 接口 |

---

## 6. 核心引擎接口（metasurface_ir 用到的最小集）

| 类型/函数 | 位置 | 说明 |
|---|---|---|
| `AbstractMaterial` + 7 子类型 | `src/Material.jl` | 统一 `𝜀̃(m, f)`（介电常数）、`ñ(m, f)`（折射率）接口 |
| `HomogeneousLayer(mat, d)` | `src/Layer.jl` | 均匀层，`d=Inf` 表半空间 |
| `PeriodicLayer(mat, d, pattern)` | `src/Layer.jl` | 周期图案层 |
| `FourierPeriodicPattern(gvs, coeffs)` | `src/Layer.jl` | 傅里叶系数形式图案，`coeffs` 为常数或 `f→Complex` 函数 |
| `RCWAGrid(b1, b2; Nx, Ny)` | `src/RCWA.jl` | 倒格子基 + 截断阶，谐波数 `nharmonics = (2Nx+1)(2Ny+1)` |
| `convolution_matrix(layer, grid, f)` | `src/RCWA.jl` | 卷积矩阵 `[ε_r]_{pq} = ε̂(g_p−g_q)` |
| `smatrix_k(layers, grid, f, kx, ky)` | `src/RCWA.jl` | 散射矩阵（内部对均匀层走 block-diag 快路径） |
| `star(SA, SB)` | `src/RCWA.jl` | Redheffer 星积 |
| `solve(layers, grid, f; θ, φ, pol_TE, pol_TM)` | `src/RCWA.jl` | RCWA 求解，返回 `RCWAResult(R, T, Rtot, Ttot, A)` |
| `solve(layers, f; θ)` | `src/TMM.jl` | 同名 TMM 重载（全均匀层时分派到这里） |

**时间约定**：`e^{-iωt}`；频率 `f` 单位 **Hz**；长度/波长单位 **米**；`ñ = n + ik`（`Im ε ≥ 0` 表示损耗）。

---

## 7. 单位约定（重要，勿混淆）

| 物理量 | 单位 |
|---|---|
| 长度 / 波长 / 厚度 | 米 (m) |
| 频率 `f` | Hz |
| Drude `omega_p` / `gamma_p` | rad/s（角频率） |
| Lorentz `f0` / `gamma` | Hz（内部 `config.jl` 自动 ×2π → rad/s） |
| Debye `tau` | 秒 (s) |
| 角度 `theta` / `phi` | 弧度 (rad) |
| `SampledDataMaterial` 波长 | **nm**，且必须升序 |
| CSV 波长列 | µm（`load_sampled_nk_csv` 默认 `λ_unit=1e-6`，转成 nm） |

---

## 8. 历史坑与修复记录（务必保留）

1. **`SampledDataMaterial` 构造 bug（已修）**：早期该类型构造函数参数个数错误、插值模式与坐标
   网格不匹配。现统一用 `Gridded(Linear())` 线性插值（对实测 n,k 最稳健，避免样条过冲）。
   注意其波长单位是 **nm 且升序**，因此 CSV（µm、降序）必须经 `load_sampled_nk_csv` 转换。

2. **`shape_image` 的归属（已修）**：它必须定义在 `metasurface_ir.jl`（核心库），因为
   `config.jl` 与 `api.jl` 都依赖它。历史上放 `api.jl` 导致 `config.jl` 单独 include
   `metasurface_ir.jl` 时 `UndefVarError: shape_image not defined`。

3. **`_as_matrix` 已矩阵化输入（已修）**：`config.jl` 的 `_as_matrix` 必须先判断
   `x isa AbstractMatrix` 直接返回，否则会把已矩阵化的 Julia Dict 输入错误 flatten 成 `(N²,1)`。

4. **juliacall 多线程（关键）**：`juliacall` **不识别** `JULIA_NUM_THREADS`，必须用
   `PYTHON_JULIACALL_THREADS`（+ `PYTHON_JULIACALL_HANDLE_SIGNALS=yes`），且在
   `import juliacall` **之前**设置。Julia 侧再 `BLAS.set_num_threads(1)` 避免与
   `Threads.@threads` 扫频 oversubscription。已在 `gestalt_metasurface.py` 顶部处理。

5. **Windows GBK 控制台**：Python 打印含 µ/中文时会因 GBK 编码崩溃。`example.py` 用
   `sys.stdout.reconfigure(encoding="utf-8", errors="replace")` 规避。

6. **Windows 下 notebook 冒烟测试**：`demo.ipynb` 不能用 `nbclient` 冒烟测试（`zmq` 事件循环
   在 Windows 会 hang）。改为自定义脚本逐个 `exec` 代码 cell + `matplotlib` 后端设为 `Agg`。

7. **block-diagonal 多层快路径（v2 优化）**：均匀层本征模态满足 `W=I`、`V=Q·λ⁻¹` 分块对角
   （每个平面波谐波一个独立 2×2 (x,y) 偏振子块），故均匀层/界面走 `BlockDiagS` 的
   `O(NH)` 星积，只有图案层走稠密 `O(NH³)`。`smatrix_k` 内部自动选择路径，接口不变。

8. **能量守恒由构造保证**：`A = 1 − R − T`，`R+T+A=1` 恒成立（约 1e-16），不要另算 A。

---

## 9. 验证与测试状态

| 脚本 | 验证内容 | 结论 |
|---|---|---|
| `validate_rcwa.jl` | 无损耗图案 `R+T=1` | max|A| ≈ 1e-12 |
| 〃 | 均匀 VO2\|Al：RCWA ↔ TMM | ΔR/ΔT ~ 1e-14 |
| 〃 | 吸收谱随 N 收敛 | N≥6 相对 N=8 < 1% |
| `test_multilayer.jl` | 新旧路径 S 矩阵逐元素对比 | ~1e-9（舍入） |
| 〃 | 纯均匀/混合/多图案层 + 性能 | block-diag 加速 ~23–37× |
| `test_config.jl` | 材料工厂 / JSON 往返 / 灰度+多层 | 通过 |
| `python/test_config.py` | `simulate` / `simulate_from_json` / `save_result` | 通过 |

---

## 10. 性能特征（benchmark.jl，16 线程）

| 参数 | 影响 | 实测 |
|---|---|---|
| `N`（截断阶） | **主导，~O(N⁶)** | 单 solve：N=4→68 ms，N=6→475 ms，N=8→2117 ms |
| `n_lambda` | 线性 | N=4 时 ~32 ms/波长 |
| 图案层数（`nslices`） | 线性 | N=6 时 ~+300 ms/片 |
| 均匀层 | 几乎免费 | block-diag，比图案层便宜 ~37× |

批量生成推荐：`N=4`（已收敛 <1%）、`n_lambda=64–100`、二值图案（避免阶梯切片）、16 线程。

---

## 11. 如何快速重新跑起来

```bash
# Julia：验证 / 形状可视化 / 配置测试 / 性能基准
julia --project=. -t auto metasurface_ir/validate_rcwa.jl
julia --project=. -t auto metasurface_ir/run_shapes.jl
julia --project=. -t auto metasurface_ir/test_config.jl
julia --project=. -t auto metasurface_ir/test_multilayer.jl
julia --project=. -t auto metasurface_ir/benchmark.jl

# Python：先装依赖
pip install juliacall numpy matplotlib jupyter
cd metasurface_ir/python && python test_config.py && python example.py
```

> 首次运行会预编译 GESTALT（数十秒）。JSON 依赖 `JSON3`（Julia 侧需 `Pkg.add("JSON3")`，
> 若报 `using JSON3` 未找到）。

---

## 12. 已知待办 / 未来方向

- Li 因式分解（高对比介电常数的收敛性，当前直接用 Fourier 系数）。
- 跨频率缓存卷积矩阵（目前每个频率重算一次 `convolution_matrix`）。
- 批量生成驱动（在样本维度继续并行）。
- 把 block-diag 快路径扩展到扭转 RCWA4D。
- 更多图案基元 + 任意晶格类型。
