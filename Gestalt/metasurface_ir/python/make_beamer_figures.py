# -*- coding: utf-8 -*-
"""Generate figures for the metasurface_ir beamer slides (work-progress report).

Figures are produced through the actual Python API so they faithfully reflect the
workflow. Benchmark scaling data is taken from `benchmark.jl` (measured on a
16-thread machine); it is embedded here to avoid re-running the benchmark.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FIGDIR = os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(FIGDIR, exist_ok=True)

from gestalt_metasurface import run_shape, shape_image, simulate

plt.rcParams.update({"font.size": 11, "figure.dpi": 110})

# ---- 1. Built-in shape patterns ------------------------------------------------
shapes = ["rectangle", "circle", "triangle", "cross"]
fig, axes = plt.subplots(1, 4, figsize=(10, 2.8))
for a, s in zip(axes, shapes):
    a.imshow(shape_image(s, N=64), cmap="viridis", origin="lower")
    a.set_title(s); a.set_xticks([]); a.set_yticks([])
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "ms_shapes.png"))
plt.close(fig)

# ---- 2. Circle absorption (2-layer stack), both phases -------------------------
res = {ph: run_shape("circle", ph, N=4, n_lambda=80) for ph in ("insulating", "metallic")}
fig, ax = plt.subplots(figsize=(5.2, 3.2))
for ph, c in (("insulating", "tab:blue"), ("metallic", "tab:red")):
    ax.plot(res[ph]["wavelength"] * 1e6, res[ph]["A"], color=c, label=ph)
ax.set(xlabel="wavelength (um)", ylabel="Absorptance A", ylim=(0, 1),
       title="circle: air | VO2(4um) | Al(200nm) | air")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "ms_absorption.png"))
plt.close(fig)

# ---- 3. API cross-check: run_shape vs simulate on the *same* stack -------------
def same_stack_config():
    return {
        "period": 5e-6, "N": 4,
        "wavelength": {"min": 8e-6, "max": 14e-6, "n": 80},
        "incidence": {"theta": 0.0, "polarization": "unpolarized"},
        "materials": {
            "air":   {"type": "nk", "n": 1.0, "k": 0.0},
            "vo2_m": {"type": "sampled_nk", "file": "data/VO2-M.csv",
                      "wavelength_unit": "micron"},
            "al":    {"type": "drude", "omega_p": 2.4e16, "gamma_p": 1.24e14,
                      "eps_inf": 1.0},
        },
        "layers": [
            {"material": "air", "thickness": None},
            {"material": "vo2_m", "pattern": {"type": "shape", "shape": "circle",
                                               "thickness": 4e-6, "background": "air"}},
            {"material": "al", "thickness": 200e-9},
            {"material": "air", "thickness": None},
        ],
    }

a = run_shape("circle", "metallic", N=4, n_lambda=80)   # convenience API
b = simulate(same_stack_config())                        # generic JSON API
wA = np.asarray(a["wavelength"]) * 1e6
AA = np.asarray(a["A"])
AB = np.asarray(b["A"])
dA = AA - AB

fig, (ax, ax2) = plt.subplots(2, 1, figsize=(5.2, 4.0),
                              gridspec_kw={"height_ratios": [2, 1]}, sharex=True)
ax.plot(wA, AA, color="tab:red", lw=2, label="run_shape (convenience)")
ax.plot(wA, AB, color="k", ls="--", lw=1.5, label="simulate (JSON config)")
ax.set(ylabel="Absorptance A", ylim=(0, 1), title="same circle stack, two APIs")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax2.plot(wA, dA, color="tab:green")
ax2.axhline(0, color="k", lw=0.8)
ax2.set(xlabel="wavelength (um)", ylabel="A difference")
ax2.text(0.02, 0.7, f"max|A_run_shape - A_simulate| = {np.max(np.abs(dA)):.2e}",
         transform=ax2.transAxes, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "ms_api_crosscheck.png"))
plt.close(fig)

# ---- 4. Benchmark scaling (from benchmark.jl, 16 threads) ----------------------
N   = np.array([2, 3, 4, 5, 6, 7, 8])
t_s = np.array([4.04, 21.28, 67.79, 198.49, 475.2, 1015.07, 2117.06])      # ms/solve
nl  = np.array([20, 50, 100, 200, 400])
t_sw= np.array([0.911, 1.788, 3.25, 6.377, 12.728])                        # s/sweep
ns  = np.array([1, 2, 4, 8, 16])
t_ns= np.array([443.6, 741.8, 1350.1, 2557.0, 4954.7])                     # ms/solve @N=6

fig, axes = plt.subplots(1, 3, figsize=(11, 3.0))
axes[0].plot(N, t_s, "o-", color="tab:red")
axes[0].set_yscale("log")
axes[0].set(xlabel="N (truncation)", ylabel="ms / solve", title="single solve vs N")
axes[0].grid(alpha=0.3, which="both")

axes[1].plot(nl, t_sw, "o-", color="tab:blue")
axes[1].set(xlabel="n_lambda", ylabel="s / sweep", title="sweep vs n_lambda (N=4)")
axes[1].grid(alpha=0.3)

axes[2].plot(ns, t_ns, "o-", color="tab:green")
axes[2].set(xlabel="nslices (# patterned layers)", ylabel="ms / solve",
            title="slices vs time (N=6)")
axes[2].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "ms_scaling.png"))
plt.close(fig)

print("Figures written to", FIGDIR)
for f in sorted(os.listdir(FIGDIR)):
    print("  ", f)
