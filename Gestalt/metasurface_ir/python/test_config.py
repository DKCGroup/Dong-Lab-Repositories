# -*- coding: utf-8 -*-
"""配置驱动通用接口的 Python 端测试。"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gestalt_metasurface import simulate, simulate_from_json, save_result

# ---- 1. dict 配置（shape 图案 + 多层） --------------------------------------
print("== 1. simulate(dict) ==")
config = {
    "period": 5e-6,
    "N": 4,
    "wavelength": {"min": 8e-6, "max": 14e-6, "n": 6},
    "incidence": {"theta": 0.0, "polarization": "unpolarized"},
    "materials": {
        "air":     {"type": "nk", "n": 1.0, "k": 0.0},
        "vo2_ins": {"type": "sampled_nk", "file": "data/VO2-I.csv", "wavelength_unit": "micron"},
        "al":      {"type": "drude", "omega_p": 2.4e16, "gamma_p": 1.24e14, "eps_inf": 1.0},
    },
    "layers": [
        {"material": "air", "thickness": None},
        {"material": "vo2_ins", "pattern": {"type": "shape", "shape": "circle",
                                             "thickness": 4e-6, "background": "air"}},
        {"material": "al", "thickness": 200e-9},
        {"material": "air", "thickness": None},
    ],
}
res = simulate(config)
print("  keys:", sorted(res.keys()))
print("  polarization:", res["polarization"], " n_harmonics:", res["n_harmonics"])
print("  A:", [round(a, 3) for a in res["A"]])

# ---- 2. JSON 文件往返 -------------------------------------------------------
print("\n== 2. simulate_from_json / save_result ==")
with open("_cfg.json", "w", encoding="utf-8") as f:
    json.dump(config, f)

res2 = simulate_from_json("_cfg.json")
save_result(res2, "_res.json")
with open("_res.json", "r", encoding="utf-8") as f:
    reloaded = json.load(f)
print("  reloaded A:", [round(a, 3) for a in reloaded["A"]])
print("  match:", reloaded["A"] == res2["A"])

os.remove("_cfg.json"); os.remove("_res.json")
print("\nALL OK")
