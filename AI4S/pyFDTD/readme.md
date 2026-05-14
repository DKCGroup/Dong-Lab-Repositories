# pyFDTD

Python helpers to build and drive **ANSYS Lumerical FDTD** simulations from code. The core module `pyFDTD.py` wraps `lumapi` with a small object-oriented layer (geometry, FDTD region, sources, monitors, meshes). The repository also includes **TARC** and **MOEMS** workflows that turn a grayscale layout image into GDS, import it into FDTD, run the solver, and return spectra or circular-dichroism–related quantities.

## Requirements

- **ANSYS Lumerical FDTD** (or compatible install) with the **Python API** (`lumapi`). Typical API path on Windows:

  `C:\Program Files\ANSYS Inc\v252\Lumerical\api\python`

  Adjust the version folder (`v252`, etc.) to match your installation.

- **Python 3.11** with at least:
  - `numpy`
  - `Pillow` (`PIL`)
  - `matplotlib` (for the example drivers)
  - `gdsfactory` (for image → GDS in `TARC.py` / `MOEMS.py`)

- **`TARC_material.mdf`** in the project root (or the working directory from which you run scripts). The TARC/MOEMS scripts copy this file into each example folder before `importmaterialdb`. This file is **not** shipped in the repository; add your own material database or obtain it from your group’s data store.

## Setup

1. Clone or copy this repository.
2. Ensure Python can import `lumapi`, for example by prepending the Lumerical API path (as in `TARC.py` / `MOEMS.py`):

   ```python
   import sys
   sys.path.append(r"C:\Program Files\ANSYS Inc\v252\Lumerical\api\python")
   ```

   Alternatively set `PYTHONPATH` to that directory in your environment.

3. Place `TARC_material.mdf` next to the scripts if you use `TARC` / `MOEMS` / the batch examples.

## Repository layout

| File | Role |
|------|------|
| `pyFDTD.py` | `Model`, `StructureGroup`, `Element`, `Sim`, `Source`, `Monitor`, `Mesh`; Palik-style material name shortcuts. |
| `TARC.py` | `TARC_generate_by_image`: MIR band (8–13 µm), absorption vs wavelength for VO₂ phases **M** / **I**. |
| `MOEMS.py` | `MOEMS_generate_by_image`: ~10.6 µm band, LCP/RCP sources for CD-style post-processing. |
| `MOEMS_generate_by_image.py` | Variant MOEMS pipeline (similar dependencies). |
| `example_run_TARC.py` | Loops over `examples/example_*`, runs M and I, saves `*_spectrum.png`. |
| `example_run_MOEMS.py` | Loops over examples, runs LCP/RCP for phase M, saves `*_CD.png`. |
| `examples/` | Per-case folders; each expects `{folder_name}.png` (and will write `.gds`, `.fsp`, logs). |

## Quick start

From the repository root, with `TARC_material.mdf` present and Lumerical licensed:

```bash
python example_run_TARC.py
```

or

```bash
python example_run_MOEMS.py
```

Each subdirectory under `examples/` should be named consistently with its input image (e.g. `examples/example_1/example_1.png`).

## Using `pyFDTD.Model` in your own script

```python
import pyFDTD

fdtd = pyFDTD.Model(r"D:\work\my_sim\", "my_case", show=False)
fdtd.switch_to_layout()
fdtd.import_material("my_materials.mdf")
fdtd.set_global_monitor("MIR")   # or "10.6um"
fdtd.set_global_source("MIR")

rect = fdtd.add_rect("substrate")
rect.set_by_min_max("x", 0, 5e-6)
rect.set_by_min_max("y", 0, 5e-6)
rect.set_by_min_max("z", -1e-6, 0)
rect.set_material("Si")
rect.assemble()

sim = fdtd.add_fdtd("fdtd")
# ... configure sim, sources, monitors ...
sim.assemble()

fdtd.save()
fdtd.run()
data = fdtd.get_result("monitor_name", "field_name")
fdtd.close()
```

Extend the classes in `pyFDTD.py` as needed; some `Source` / `Monitor` branches are marked TBD for advanced options.