# 红外超表面（IR Metasurface）工作流

面向 **8–14 µm 长波红外（LWIR）** 的超表面仿真工作流，基于 GESTALT 的 RCWA 求解器，
同时提供 **Julia** 与 **Python**（juliacall）两套调用接口。

## 目录结构

```
metasurface_ir/
├── metasurface_ir.jl        # 核心库：材料 / 灰度图→图案 / 堆栈构建 / 扫频 / 形状生成
├── config.jl                # 配置驱动通用求解器（JSON schema → simulate）
├── api.jl                   # 跨语言 API（Python 调用入口，ASCII 参数名，返回纯数据）
├── validate_rcwa.jl         # RCWA 正确性验证（能量守恒 / RCWA↔TMM / 收敛性）
├── test_config.jl           # 配置驱动接口回归测试（材料/图案/JSON 往返）
├── test_multilayer.jl       # 多层快速路径正确性 + 性能基准
├── benchmark.jl             # 批量生成性能基准（N / n_lambda / nslices 耗时）
├── run_shapes.jl            # 形状仿真 + 可视化（矩形/圆/三角/十字）
├── figures/                 # 输出图像
├── README.md                # 使用者文档（JSON schema / API / 性能）
├── ARCHIVE.md               # 维护者归档（代码结构 / 数据流 / 历史坑 / 验证）
└── python/
    ├── gestalt_metasurface.py   # Python API 封装（含 simulate 通用接口）
    ├── example.py               # Python 使用示例
    ├── demo.ipynb               # Jupyter 演示 notebook（展示完整运行流程）
    ├── test_config.py           # Python 端配置接口回归测试
    └── requirements.txt         # Python 依赖

data/                        # 实测材料数据（项目根目录）
├── VO2-I.csv                # 绝缘态（monoclinic）VO2 的 n,k
└── VO2-M.csv                # 金属态（rutile）VO2 的 n,k
```

## 结构定义

工作波段 8–14 µm，周期 `p = 5 µm` 的正方晶格，标准两层结构：

```
入射：空气（半无限）
顶层：VO2 图案层，厚度由灰度图决定（0–255 → 0–4 µm）
底层：铝（Al）衬底，厚度固定 200 nm
出射：空气（半无限）
```

- **灰度图输入**：`N×N`（默认 64×64）灰度矩阵，取值 0–255。
  - 灰度 0（或 ≤ 阈值 `t_min_gray`，默认 50）→ 空气；
  - 灰度 `g > 50` → VO2，厚度 = `(g/255) × max_t`（`max_t` 默认 4 µm）。
  - **二值图**（0 与某个 >50 的值）→ 单层 `PeriodicLayer`（RCWA 快路径）；
  - **多高度灰度图** → 自动按 `nslices`（默认 8）个等厚子层做阶梯（staircase）切片。
- **多层扩展预留**：`build_stack(..., spacers=...)` 可在图案层与铝镜之间插入任意数量的
  均匀无图案介电层（`HomogeneousLayer`）。
- **底部铝镜**：200 nm 铝在 LWIR 约十余个趋肤深度，**完全不透明**（T ≈ 0），因此这是
  一个 **吸收器/反射器**，A = 1 − R。若要非零透射，需减薄/去除铝镜或改用透明衬底。

## VO2 材料数据

`data/VO2-I.csv` 与 `data/VO2-M.csv` 是实测 VO2 折射率表，三列：

```
wavelengthmicrons, Reindex, Imindex     # 波长(µm), n, k
```

覆盖 0.4–30 µm（数据为降序排列）。接入时由 `load_sampled_nk_csv` 自动完成：
**降序→升序**、**微米→纳米**（GESTALT 的 `SampledDataMaterial` 波长单位为 nm）、
**去重**。默认采用 **线性插值**（对实测 n,k 最稳健，避免样条过冲）。

> 注意：本工作流修复了 `SampledDataMaterial` 的一个历史 bug（构造函数参数个数错误、
> 插值模式与坐标网格不匹配），现已可用实测数据正常构造。

材料接口：

```julia
vo2_insulating_sampled()   # 实测绝缘态 VO2（data/VO2-I.csv）
vo2_metallic_sampled()     # 实测金属态 VO2（data/VO2-M.csv）
# 拟合模型（无实测数据时的回退）：
vo2_insulating_material()  # Lorentz 单振子
vo2_metallic_material()    # Drude
```

实测数据特征（用于核对）：

| 相态 | n @ 8 µm | k @ 8 µm | n @ 14 µm | k @ 14 µm |
|---|---|---|---|---|
| 绝缘态（I） | 2.53 | 0.048 | 0.46 | 1.69 |
| 金属态（M） | 6.22 | 9.48 | 9.50 | 12.57 |

---

# 一、Julia 用法

## 运行验证与示例

```bash
# RCWA 正确性验证（能量守恒 / RCWA↔TMM 交叉 / 收敛性）
julia --project=. -t auto metasurface_ir/validate_rcwa.jl

# 形状仿真 + 可视化
julia --project=. -t auto metasurface_ir/run_shapes.jl
```

## 最小示例

```julia
using GESTALT
include("metasurface_ir/metasurface_ir.jl")
include("metasurface_ir/api.jl")

# 圆形图案，金属态 VO2，8–14 µm 光谱
res = metasurface_shape_spectrum("circle", "metallic"; n_lambda=200, N=6)
res.wavelength   # 波长（米）
res.A            # 吸收谱
```

---

# 二、Python 用法（合作者）

通过 [juliacall](https://github.com/JuliaPy/PythonCall.jl) 调用 GESTALT，Python 为主调方，
Julia 为计算引擎。

## 环境准备

1. 安装 **Julia**（系统 PATH 中需有 `julia`）。
2. 安装 Python 依赖：

```bash
pip install juliacall numpy matplotlib jupyter
```

> **Jupyter 演示**：`python/demo.ipynb` 是一个可直接运行的演示 notebook，
> 逐步展示色散、形状图案、光谱与通用多层接口。在 `metasurface_ir/python/` 目录下执行
> `jupyter notebook`（或 `jupyter lab`）后打开即可。

## 快速开始

将 `metasurface_ir/python/` 加入 `sys.path`（或把 `gestalt_metasurface.py` 复制到你的项目）：

```python
import sys
sys.path.insert(0, r"C:\Research\GESTALT\metasurface_ir\python")

import numpy as np
from gestalt_metasurface import run_shape, run_spectrum, vo2_dispersion
from gestalt_metasurface import simulate, simulate_from_json, save_result

# 1) 内置形状 → R/T/A 光谱
res = run_shape("circle", "metallic", n_lambda=256, N=6)
res["wavelength"]   # numpy 数组，单位米
res["R"], res["T"], res["A"]

# 2) 自定义灰度图（64×64，0 = 空气，>50 = VO2）
img = np.zeros((64, 64), dtype=np.uint8)
img[20:44, 20:44] = 255
res = run_spectrum(img, "insulating", N=6, n_lambda=200)

# 3) VO2 复折射率色散
d = vo2_dispersion("metallic", n_points=200)
d["wavelength"], d["n"], d["k"]

# 4) 通用配置驱动（任意多层 / 任意材料，JSON 进出）
config = {
    "period": 5e-6, "N": 6,
    "wavelength": {"min": 8e-6, "max": 14e-6, "n": 200},
    "materials": {
        "air": {"type": "nk", "n": 1.0, "k": 0.0},
        "vo2": {"type": "sampled_nk", "file": "data/VO2-M.csv", "wavelength_unit": "micron"},
        "al":  {"type": "drude", "omega_p": 2.4e16, "gamma_p": 1.24e14, "eps_inf": 1.0},
    },
    "layers": [
        {"material": "air", "thickness": None},
        {"material": "vo2", "pattern": {"type": "shape", "shape": "circle",
                                         "thickness": 4e-6, "background": "air"}},
        {"material": "al", "thickness": 200e-9},
        {"material": "air", "thickness": None},
    ],
}
res = simulate(config)                 # dict 进，dict 出
save_result(res, "result.json")        # 保存结果
```

首次调用会启动 Julia 进程并预编译 GESTALT（数十秒），之后复用同一进程。

> **多线程**：光谱扫描在 Julia 端多线程并行，线程数由 `PYTHON_JULIACALL_THREADS`
> 控制（`juliacall` 不识别 `JULIA_NUM_THREADS`），默认取 `min(CPU 核数, 8)`。
> 可在 `import gestalt_metasurface` **之前**设置，例如
> `os.environ["PYTHON_JULIACALL_THREADS"] = "16"`。

---

# 三、API 参考（Python）

## `run_spectrum(img, phase, ...)`

由灰度图计算 8–14 µm 的 R/T/A 光谱。

| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `img` | 二维数组 | — | `Np×Np` 灰度图，取值 0–255 |
| `phase` | str | `"insulating"` | `"insulating"` / `"metallic"` |
| `p` | float | `5e-6` | 周期（米） |
| `max_t` | float | `4e-6` | 灰度 255 对应最大厚度（米） |
| `t_min_gray` | float | `50` | 灰度阈值，≤ 此值视为空气 |
| `N` | int | `6` | RCWA 截断阶（谐波数 `(2N+1)²`） |
| `lambda_min` | float | `8e-6` | 波长下限（米） |
| `lambda_max` | float | `14e-6` | 波长上限（米） |
| `n_lambda` | int | `200` | 频率采样点数（≥2） |
| `theta` | float | `0.0` | 入射角（弧度） |
| `sampled` | bool | `True` | `True` 用实测 n,k；`False` 用拟合模型 |
| `d_al` | float | `200e-9` | 铝镜厚度（米） |

返回 `dict`，键 `"wavelength"`, `"R"`, `"T"`, `"A"`（numpy 数组）。

## `run_shape(shape, phase, **kwargs)`

内置形状 + 相态 → R/T/A 光谱。`shape` 为 `"rectangle"` / `"circle"` /
`"triangle"` / `"cross"`；`**kwargs` 透传给 `run_spectrum`。

## `shape_image(shape, N=64, gray=255)`

返回内置形状的灰度图（numpy 数组，`uint8`）。

## `vo2_dispersion(phase, lambda_min=8e-6, lambda_max=14e-6, n_points=200, sampled=True)`

返回 VO2 复折射率色散，`dict` 键 `"wavelength"`, `"n"`, `"k"`。

## 通用接口 `simulate(config)` / `simulate_from_json(path)` / `save_result(result, path)`

配置驱动的通用仿真接口，输入/输出均为 JSON（dict），详见
[「四、通用配置驱动接口」](#四通用配置驱动接口json-schema)。

```python
from gestalt_metasurface import simulate, simulate_from_json, save_result

res = simulate(config_dict)              # dict 配置 → dict 结果
res = simulate_from_json("cfg.json")     # JSON 文件 → dict 结果
save_result(res, "result.json")          # dict 结果 → JSON 文件
```

---

# 四、通用配置驱动接口（JSON schema）

这是**泛化后的核心接口**：任意多层、任意材料、任意波段的结构化仿真，输入与输出都是
格式化的 JSON（Python 端为等价的 dict）。旧的 `run_spectrum` / `run_shape` 仍保留，作为
「VO2 + 铝镜」这一常见场景的便捷封装。

## 4.1 单位约定（全部 SI）

| 物理量 | 单位 |
|---|---|
| 长度 / 波长 / 厚度 | 米 (m) |
| 频率 | 赫兹 (Hz) |
| Drude `omega_p` / `gamma_p` | 弧度/秒 (rad/s) |
| Lorentz `f0` / `gamma` | 赫兹 (Hz)，内部自动 ×2π |
| Debye `tau` | 秒 (s) |
| 角度 `theta` / `phi` | 弧度 (rad) |

## 4.2 输入配置（顶层字段）

| 字段 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `period` | float | `5e-6` | 正方晶格周期（米） |
| `N` | int | `6` | RCWA 截断阶，谐波数 `(2N+1)²` |
| `wavelength` | 对象或数组 | 8–14 µm | `{"min","max","n"}` 或显式数组 |
| `incidence` | 对象 | 正入射非偏振 | `{"theta","phi","polarization"}` |
| `materials` | 对象 | — | 名称 → 材料定义（见 4.3） |
| `layers` | 数组 | — | 层列表（见 4.4） |

`incidence.polarization` 取值 `"TE"` / `"TM"` / `"unpolarized"`（默认）。

## 4.3 材料定义（`materials` 内）

| `type` | 字段 | 说明 |
|---|---|---|
| `nk` | `n`, `k` | 常数折射率 / 消光系数 |
| `drude` | `omega_p`, `gamma_p`, `eps_inf` | Drude 金属（rad/s） |
| `lorentz` | `eps_inf`, `oscillators` | 振子数组 `[{"deps","f0","gamma"}]`（Hz） |
| `debye` | `eps_s`, `eps_inf`, `tau` | Debye 弛豫（秒） |
| `sampled_nk` | `file` + `wavelength_unit`，或内联 `wavelength_nm`/`n`/`k` | 实测 n,k 表 |

`sampled_nk` 两种写法：

```json
{"type": "sampled_nk", "file": "data/VO2-M.csv", "wavelength_unit": "micron"}
{"type": "sampled_nk", "wavelength_nm": [...], "n": [...], "k": [...]}
```

## 4.4 层定义（`layers` 数组）

| 字段 | 类型 | 说明 |
|---|---|---|
| `material` | str | 引用 `materials` 中的名称 |
| `thickness` | float 或 `null` | 均匀层厚度（米）；`null` = 半空间（仅首/末层） |
| `pattern` | 对象 | 可选；有 `pattern` 即图案层（厚度由图案决定） |

图案层有三种 `pattern.type`：

| `type` | 字段 | 说明 |
|---|---|---|
| `grayscale` | `image`, `max_thickness`, `min_gray`, `max_gray`, `background`, `nslices` | 灰度→厚度，多高度自动阶梯切片 |
| `binary` | `image`, `thickness`, `background` | 二值图（0/非0）+ 固定厚度 |
| `shape` | `shape`, `thickness`, `background`, `pixels`, `gray` | 内置形状（rectangle/circle/triangle/cross） |

`background` 是图案空隙处的材料名（默认 `"air"`）。首层与末层必须是半空间
（入射 / 出射介质）。

## 4.5 完整示例

```json
{
  "period": 5e-6,
  "N": 6,
  "wavelength": {"min": 8e-6, "max": 14e-6, "n": 200},
  "incidence": {"theta": 0.0, "polarization": "unpolarized"},
  "materials": {
    "air":     {"type": "nk", "n": 1.0, "k": 0.0},
    "vo2_ins": {"type": "sampled_nk", "file": "data/VO2-I.csv", "wavelength_unit": "micron"},
    "vo2_met": {"type": "sampled_nk", "file": "data/VO2-M.csv", "wavelength_unit": "micron"},
    "al":      {"type": "drude", "omega_p": 2.4e16, "gamma_p": 1.24e14, "eps_inf": 1.0},
    "ge":      {"type": "nk", "n": 4.0, "k": 0.0}
  },
  "layers": [
    {"material": "air", "thickness": null},
    {"material": "vo2_ins", "pattern": {
        "type": "grayscale", "image": [[0,0,255],[0,255,255],[255,255,0]],
        "max_thickness": 4e-6, "min_gray": 50, "max_gray": 255,
        "background": "air", "nslices": 8}},
    {"material": "ge", "thickness": 500e-9},
    {"material": "al", "thickness": 200e-9},
    {"material": "air", "thickness": null}
  ]
}
```

## 4.6 输出结果

`simulate` 返回的 dict（可直接 `json.dump`）字段：

| 键 | 类型 | 说明 |
|---|---|---|
| `wavelength` / `wavelength_um` | 数组 | 波长（米 / 微米） |
| `R` / `T` / `A` | 数组 | 组合偏振的反射 / 透射 / 吸收 |
| `R_TE` / `T_TE` / `A_TE` | 数组 | TE 偏振分量 |
| `R_TM` / `T_TM` / `A_TM` | 数组 | TM 偏振分量 |
| `polarization` | str | 请求的偏振（`unpolarized` 时 `R/T/A` 为平均） |
| `period` | float | 周期（米） |
| `n_harmonics` | int | 平面波谐波数 `(2N+1)²` |
| `n_layers` | int | 实际层数（含阶梯切片展开） |
| `max_abs_RTA_deviation` | float | 能量守恒偏差 `max|R+T+A−1|` |

> 无论请求哪种偏振，TE 与 TM 分量都会同时返回，结构一致，便于后续分析。

---

# 五、Julia API

## `api.jl`（便捷封装，参数名与 Python 一致）

```julia
metasurface_spectrum(img, phase; p, max_t, t_min_gray, N, lambda_min, lambda_max,
                     n_lambda, theta, sampled, d_al)      # → (wavelength, R, T, A)
metasurface_shape_spectrum(shape, phase; kwargs...)       # → (wavelength, R, T, A)
shape_image(shape; N, gray)                              # → Matrix{UInt8}
vo2_dispersion(phase; lambdas, sampled)                   # → (wavelength, n, k)
```

## `config.jl`（通用接口）

```julia
include("metasurface_ir/config.jl")

simulate(config)               # Dict/JSON3.Object 配置 → NamedTuple 结果
simulate_json(json_str)        # JSON 字符串 → JSON 字符串（跨语言主入口）
material_from_config(cfg)      # 材料定义 → AbstractMaterial
stack_from_config(config)      # 配置 → (grid, layers)
read_config_json(path)         # JSON 文件 → 解析对象
```

---

# 六、性能与精度

- **截断阶 N**：谐波数 `(2N+1)²`。收敛性已验证（`validate_rcwa.jl`）：
  N ≥ 6 相对 N = 8 的误差 < 1%。单次 `solve` 耗时（TE 或 TM，`benchmark.jl` 实测，
  16 线程机器）基本按 O(N⁶) 增长：
  N=4 ≈ 68 ms、N=6 ≈ 475 ms、N=8 ≈ 2117 ms。
- **能量守恒**：无损耗介质 R+T=1 满足到 ~1e-12；RCWA 与独立 TMM 求解器一致性 ~1e-14。

## 批量生成性能基准（benchmark.jl）

面向自动化批量训练集生成，`metasurface_ir/benchmark.jl` 量化了各关键参数对单次
运行时间的影响（`julia --project=. -t auto metasurface_ir/benchmark.jl`）：

| 参数 | 影响 | 实测（16 线程） |
|---|---|---|
| `N`（截断阶） | **主导**，~O(N⁶) | 单 solve：N=4→68 ms，N=6→475 ms（7×），N=8→2117 ms（31×） |
| `n_lambda`（采样点） | 线性 | N=4 时每波长 ≈ 32 ms |
| 图案层数（nslices） | 线性 | N=6 时每多一片 ≈ +300 ms |
| 均匀层 | 几乎免费 | block-diagonal 快路径，比图案层便宜 ≈ 37× |

批量生成的推荐配置（速度/精度平衡）：`N=4`（相对 N=8 已收敛 <1%）、
`n_lambda=64–100`（谱形平滑时）、二值图案（避免多高度灰度阶梯切片）、线程数 = CPU 核数。
据此估算，N=4、`n_lambda=100` 时单样本 ≈ 3.3 s；1 万样本数据集 ≈ 6–9 小时
（可在样本维度继续并行）。

## 多层快速路径（block-diagonal，v2 优化）

均匀层在 RCWA 中原本按稠密 `2NH×2NH` 矩阵逐层求逆（`O(NH³)`）。利用均匀介质的
本征模态 `W = I`、`V = Q·λ⁻¹` 是**分块对角**（每个平面波谐波一个独立的 2×2 (x,y) 偏振
子块）这一事实，`smatrix_k` 现在把**所有均匀层与均匀界面**改走 2×2 分块对角快速路径
（`O(NH)`），仅图案层保留稠密 `O(NH³)` 计算。

- 效果（`test_multilayer.jl`）：纯均匀多层 **~23–29× 加速**；图案层 + 多均匀层的
  混合结构，均匀段几乎零开销。
- 数值一致性：新旧路径 S 矩阵逐元素差 ~1e-9（双精度舍入），R/T/A 一致到 1e-12，
  能量守恒保持 `R+T+A=1`。
- 覆盖场景：任意多层均匀层（spacer / 介电 / 金属 / 衬底）、单个或多个图案层
  与均匀层任意交错，均已验证。
