# -*- coding: utf-8 -*-
"""GESTALT IR-metasurface Python API usage example.

Install deps:  pip install juliacall numpy matplotlib
Run:           python example.py

Note: with n_lambda=200 and N=6, each shape×phase sweep takes ~1–2 min
(longer on first run, which precompiles GESTALT).
"""

import sys

# 统一 stdout 为 UTF-8，避免 Windows GBK 控制台打印 Unicode（µ、中文）时崩溃。
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np

from gestalt_metasurface import run_spectrum, run_shape, shape_image, vo2_dispersion

# ---- 1. VO2 measured dispersion ---------------------------------------------
print("== 1. VO2 dispersion ==")
for phase in ("insulating", "metallic"):
    d = vo2_dispersion(phase, n_points=5)
    print(f"  {phase:>10s}  n @ 8um={d['n'][0]:.3f}  k @ 8um={d['k'][0]:.3f}")

# ---- 2/3. built-in shape spectra --------------------------------------------
print("\n== 2. circle / metallic ==")
res_met = run_shape("circle", "metallic", n_lambda=200, N=6)
print(f"  A(10um)={res_met['A'][np.argmin(abs(res_met['wavelength']-10e-6))]:.3f}")

print("\n== 3. circle / insulating ==")
res_ins = run_shape("circle", "insulating", n_lambda=200, N=6)
print(f"  A(10um)={res_ins['A'][np.argmin(abs(res_ins['wavelength']-10e-6))]:.3f}")

# ---- 4. custom grayscale image (64x64) --------------------------------------
print("\n== 4. custom grayscale image ==")
img = np.zeros((64, 64), dtype=np.uint8)
img[16:48, 24:40] = 255          # central rectangular VO2 pillar
res_custom = run_spectrum(img, "insulating", N=6, n_lambda=200)
print(f"  R range {res_custom['R'].min():.3f}-{res_custom['R'].max():.3f}, "
      f"A range {res_custom['A'].min():.3f}-{res_custom['A'].max():.3f}")

# ---- 5. visualization (optional, needs matplotlib) --------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    im = shape_image("circle", N=64)
    axes[0].imshow(im, cmap="viridis", origin="lower")
    axes[0].set_title("pattern (circle)")

    axes[1].plot(res_ins["wavelength"] * 1e6, res_ins["A"], "-", label="insulating")
    axes[1].plot(res_met["wavelength"] * 1e6, res_met["A"], "--", label="metallic")
    axes[1].set_xlabel("wavelength (um)")
    axes[1].set_ylabel("absorptance A")
    axes[1].set_title("8-14 um absorption")
    axes[1].legend()
    axes[1].set_ylim(0, 1)

    fig.tight_layout()
    fig.savefig("metasurface_python_demo.png", dpi=120)
    print("\n== 5. saved metasurface_python_demo.png ==")
except ImportError:
    print("\n== 5. matplotlib not installed, skip plotting ==")

print("\nDone.")
