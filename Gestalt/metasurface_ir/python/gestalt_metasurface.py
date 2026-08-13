"""
GESTALT — 红外超表面 Python API。

通过 juliacall（Julia↔Python 官方桥接）调用 GESTALT 的 RCWA 求解器，为 Python
用户提供红外波段（8–14 µm）超表面的仿真接口。

依赖
----
- Python 3.9+，已安装 ``juliacall``（``pip install juliacall``）
- 系统已安装 Julia（本机 PATH 中需有 ``julia``）

首次调用会启动一个 Julia 进程并预编译 GESTALT（约几秒到几十秒），之后即复用。

快速开始
--------
>>> import numpy as np
>>> from gestalt_metasurface import run_spectrum, run_shape
>>> # 用形状名直接算
>>> res = run_shape("circle", "metallic", n_lambda=256)
>>> res["A"]                     # 吸收谱 numpy.ndarray
>>> # 或用自定义灰度图（64×64，0 = 空气，>50 = VO2）
>>> img = np.zeros((64, 64), dtype=np.uint8)
>>> img[20:44, 20:44] = 255
>>> res = run_spectrum(img, "insulating", N=6)
"""

import os
import json

# 在导入 juliacall 之前设置 Julia 多线程。注意：juliacall 用 PYTHON_JULIACALL_THREADS
# 控制线程数（而非 JULIA_NUM_THREADS，后者对 juliacall 不生效）。GESTALT 的光谱扫描
# 在 Julia 端多线程并行，线程数越大扫频越快。
# HANDLE_SIGNALS 是多线程模式下 Julia GC 所必需（见 PythonCall 文档），避免段错误。
os.environ.setdefault("PYTHON_JULIACALL_THREADS", str(min(os.cpu_count() or 4, 8)))
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")

import numpy as np
import juliacall
from juliacall import Main as jl

# GESTALT 项目根目录：本文件位于 <root>/metasurface_ir/python/ 下。
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_METASURFACE_DIR = os.path.join(_ROOT, "metasurface_ir")

_loaded = False


def _ensure_loaded():
    """加载 GESTALT 与 IR 超表面工作流（幂等）。"""
    global _loaded
    if _loaded:
        return
    # Julia 字符串字面量用双引号；Windows 反斜杠转正斜杠避免转义问题。
    root = _ROOT.replace("\\", "/")
    jl.seval(f'push!(LOAD_PATH, "{root}")')
    jl.seval("using GESTALT")
    # 避免 BLAS 多线程与 Julia 多线程 oversubscription（sweep 已用 Threads.@threads 并行）
    jl.seval("using LinearAlgebra; BLAS.set_num_threads(1)")
    jl.include(os.path.join(_METASURFACE_DIR, "metasurface_ir.jl"))
    jl.include(os.path.join(_METASURFACE_DIR, "api.jl"))
    jl.include(os.path.join(_METASURFACE_DIR, "config.jl"))
    _loaded = True


def _pack(res):
    """把 Julia 返回的 NamedTuple 转成 Python dict（值为 numpy 数组）。"""
    return {
        "wavelength": np.asarray(res.wavelength, dtype=np.float64),
        "R": np.asarray(res.R, dtype=np.float64),
        "T": np.asarray(res.T, dtype=np.float64),
        "A": np.asarray(res.A, dtype=np.float64),
    }


def run_spectrum(img, phase="insulating", *, p=5e-6, max_t=4e-6, t_min_gray=50.0,
                 N=6, lambda_min=8e-6, lambda_max=14e-6, n_lambda=200,
                 theta=0.0, sampled=True, d_al=200e-9):
    """由灰度图像计算 8–14 µm 的 R/T/A 光谱。

    参数
    ----
    img : array_like, 二维 (Np×Np)
        灰度图案，取值 0–255。0 表示空气；大于 ``t_min_gray`` 的像素为 VO2，
        其厚度 = (灰度/255) × ``max_t``。二值图（0 与某个 >50 的值）即单层图案。
    phase : str
        VO2 相态，``"insulating"``（绝缘）或 ``"metallic"``（金属）。
    p : float
        超表面周期（米），默认 5 µm。
    max_t : float
        灰度 255 对应的图案层最大厚度（米），默认 4 µm。
    t_min_gray : float
        灰度阈值，≤ 该值视为空气，默认 50。
    N : int
        RCWA 平面波截断阶数（总谐波数 (2N+1)²）。N 越大越精确但越慢，默认 6。
    lambda_min / lambda_max : float
        波长范围（米），默认 8–14 µm。
    n_lambda : int
        频率采样点数（≥2），默认 200。
    theta : float
        入射角（弧度），默认 0（正入射）。
    sampled : bool
        True 使用实测 n,k 表；False 用 Lorentz/Drude 拟合模型。
    d_al : float
        底部铝镜厚度（米），默认 200 nm（在 LWIR 近似不透明）。

    返回
    ----
    dict，键为 "wavelength"、"R"、"T"、"A"（均为 numpy.ndarray，长度 n_lambda）。
    """
    _ensure_loaded()
    imgf = np.asarray(img, dtype=np.float64)
    if imgf.ndim != 2 or imgf.shape[0] != imgf.shape[1]:
        raise ValueError(f"img 必须是方阵（N×N），实际 {imgf.shape}")
    res = jl.metasurface_spectrum(
        imgf, phase,
        p=float(p), max_t=float(max_t), t_min_gray=float(t_min_gray),
        N=int(N), lambda_min=float(lambda_min), lambda_max=float(lambda_max),
        n_lambda=int(n_lambda), theta=float(theta), sampled=bool(sampled),
        d_al=float(d_al),
    )
    return _pack(res)


def shape_image(shape, N=64, gray=255):
    """生成内置形状的灰度图（0 / gray 二值）。

    参数
    ----
    shape : str
        ``"rectangle"`` / ``"circle"`` / ``"triangle"`` / ``"cross"``。
    N : int
        像素分辨率（N×N），默认 64。
    gray : int
        图案区域的灰度值（0–255），默认 255。

    返回
    ----
    numpy.ndarray，形状 (N, N)，dtype uint8。
    """
    _ensure_loaded()
    return np.asarray(jl.shape_image(str(shape), N=int(N), gray=int(gray)))


def run_shape(shape, phase="insulating", **kwargs):
    """内置形状 + 相态 → R/T/A 光谱（参数同 :func:`run_spectrum`）。

    参数
    ----
    shape : str
        ``"rectangle"`` / ``"circle"`` / ``"triangle"`` / ``"cross"``。
    phase : str
        ``"insulating"`` 或 ``"metallic"``。
    **kwargs
        透传给 :func:`run_spectrum`（p、N、n_lambda 等）。

    返回
    ----
    dict，键为 "wavelength"、"R"、"T"、"A"。
    """
    _ensure_loaded()
    res = jl.metasurface_shape_spectrum(str(shape), phase, **kwargs)
    return _pack(res)


def vo2_dispersion(phase="insulating", *, lambda_min=8e-6, lambda_max=14e-6,
                   n_points=200, sampled=True):
    """返回 VO2 的复折射率 n + i·k 随波长的色散曲线。

    返回
    ----
    dict，键为 "wavelength"、"n"、"k"（numpy.ndarray）。
    """
    _ensure_loaded()
    lambdas = np.linspace(float(lambda_min), float(lambda_max), int(n_points))
    res = jl.vo2_dispersion(phase, lambdas=lambdas, sampled=bool(sampled))
    return {
        "wavelength": np.asarray(res.wavelength, dtype=np.float64),
        "n": np.asarray(res.n, dtype=np.float64),
        "k": np.asarray(res.k, dtype=np.float64),
    }


# ---- 配置驱动（通用）接口 ---------------------------------------------------------------
# 输入是一个 JSON/dict 配置，输出是一个结构化 dict（可直接 json.dump）。单位全部 SI：
# 长度/波长/厚度 = 米，频率 = Hz，角频率 ωₚ/γ = rad/s，Lorentz f₀/γ = Hz，角度 = 弧度。

def simulate(config):
    """由配置 dict 运行一次仿真，返回结果 dict。

    参数
    ----
    config : dict
        仿真配置（JSON schema 见 README「通用配置驱动接口」）。顶层字段：
        ``period``、``N``、``wavelength``（``{"min","max","n"}`` 或显式数组）、
        ``incidence``（``{"theta","phi","polarization"}``）、
        ``materials``（名称→材料定义）、``layers``（层列表）。

    返回
    ----
    dict，键包含 "wavelength"（米）、"wavelength_um"、"R"/"T"/"A"（组合偏振）、
    "R_TE"/"T_TE"/"A_TE"、"R_TM"/"T_TM"/"A_TM"、"polarization"、
    "n_harmonics"、"n_layers"、"max_abs_RTA_deviation"。
    """
    _ensure_loaded()
    return json.loads(jl.simulate_json(json.dumps(config)))


def simulate_from_json(path):
    """从 JSON 配置文件读取并运行仿真，返回结果 dict。

    参数
    ----
    path : str
        JSON 配置文件路径。

    返回
    ----
    dict，同 :func:`simulate`。
    """
    with open(path, "r", encoding="utf-8") as f:
        config = json.load(f)
    return simulate(config)


def save_result(result, path):
    """把仿真结果 dict 保存为 JSON 文件（UTF-8，缩进 2）。"""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    return path
