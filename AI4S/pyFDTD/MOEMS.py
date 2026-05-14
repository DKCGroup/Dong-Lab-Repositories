import sys
import os
import shutil
sys.path.append(r"C:\Program Files\ANSYS Inc\v252\Lumerical\api\python")
import lumapi
import numpy as np
from PIL import Image
import gdsfactory as gf
gf.gpdk.PDK.activate()
import pyFDTD

def MOEMS_generate_by_image(filepath, filename, source="LCP", phase="M", period=5e-6, h_Pt=50e-9, h_VO2=0.16e-6, h_Cr=30e-9, h_Au=90e-9, h_air=1e-6):

    script_dir = os.path.dirname(os.path.abspath("."))
    target_dir = os.path.join(script_dir, filepath)
    os.makedirs(target_dir, exist_ok=True)
    save_path = os.path.join(target_dir, f"{filename}.fsp")
    shutil.copy("./TARC_material.mdf", f"{filepath}TARC_material.mdf")

    # Open image file
    img = Image.open(filepath+filename+".png").convert('L')
    img_arr = np.array(img)
    pixels = np.array([[0 if img_arr[i,j] < np.mean(img_arr) else 255 for j in range(img_arr.shape[1])] for i in range(img_arr.shape[0])])
    img_x, img_y = img_arr.shape
    pixel_size_x = period/img_x
    pixel_size_y = period/img_y

    # Change image to gds
    c = gf.Component()
    for i in range(img_x):
        for j in range(img_y):
            if pixels[i,j] == 0:
                c1 = c << gf.components.rectangle(size=(np.ceil(pixel_size_x*1e9)/1e3, np.ceil(pixel_size_y*1e9)/1e3), layer=(1, 0))
                c1.move((i*pixel_size_x*1e6, j*pixel_size_y*1e6))

    c.flatten()
    c.write_gds(filepath+filename+".gds")

    # Start up FDTD Model
    fdtd = pyFDTD.Model(filepath, filename+f"_{source}"+f"_{phase}")
    fdtd.switch_to_layout()
    fdtd.import_material("TARC_material.mdf")
    fdtd.set_global_monitor("10.6um")
    fdtd.set_global_source("10.6um")

    # add Pt
    Pt = fdtd.add_rect("Pt")
    Pt.set_by_min_max('x', 0, period)
    Pt.set_by_min_max('y', 0, period)
    Pt.set_by_min_max('z', -h_air-h_Pt, -h_air)
    Pt.set_material("Pt")
    Pt.set_mesh_order(1)
    Pt.assemble()

    # add Cantilever

    VO2 = fdtd.import_gds("VO2", filename, c.name, f"{c.layers[0][0]}/{c.layers[0][1]}", f"VO2-{phase}", 0, h_VO2)
    Cr = fdtd.import_gds("Cr", filename, c.name, f"{c.layers[0][0]}/{c.layers[0][1]}", "Cr", h_VO2, h_VO2+h_Cr)
    Au = fdtd.import_gds("Au", filename, c.name, f"{c.layers[0][0]}/{c.layers[0][1]}", "Au", h_VO2+h_Cr, h_VO2+h_Cr+h_Au)

    # add FDTD
    sim = fdtd.add_fdtd("fdtd")
    sim.set_by_min_max('x', 0, period)
    sim.set_by_min_max('y', 0, period)
    sim.set_by_min_max('z', -5e-6, 30e-6)
    sim.set_boundary_condition("x", "Period")
    sim.set_boundary_condition("y", "Period")
    sim.set_boundary_condition("z", "PML")
    sim.set_simulation_time(10000)
    sim.set_auto_shutoff(1e-5)
    sim.set_mesh_type("custom")
    sim.set_mesh_definition("x", "mesh step")
    sim.set_mesh_definition("y", "mesh step")
    sim.set_mesh_definition("z", "mesh step")
    sim.set_maximum_mesh_step("x", 0.1)
    sim.set_maximum_mesh_step("y", 0.1)
    sim.set_maximum_mesh_step("z", 0.1)
    sim.set_mesh_cell_per_wavelength(10)
    sim.set_dt_stability_factor(0.95)
    sim.set_minimum_mesh_step(0.001)
    sim.assemble()

    # add mesh
    mesh = fdtd.add_mesh("mesh")
    mesh.set_by_min_max("x", 0, period)
    mesh.set_by_min_max("y", 0, period)
    mesh.set_by_min_max("z", 0, h_VO2+h_Cr+h_Au)
    mesh.set_override_direction("z")
    mesh.set_override_mesh_step("z", 5e-3)
    mesh.assemble()

    # add source
    if source == "LCP":
        inc1 = fdtd.add_source("Inc1", "plane")
        inc1.set_axis_direction("-z")
        inc1.set_phase(0)
        inc1.set_by_min_max('x', 0, period)
        inc1.set_by_min_max('y', 0, period)
        inc1.set_by_value('z',25e-6)
        inc1.assemble()
        inc2 = fdtd.add_source("Inc2", "plane")
        inc2.set_axis_direction("-z")
        inc2.set_phase(90)
        inc2.set_polarization(90)
        inc2.set_by_min_max('x', 0, period)
        inc2.set_by_min_max('y', 0, period)
        inc2.set_by_value('z',25e-6)
        inc2.assemble()

    elif source == "RCP":
        inc1 = fdtd.add_source("Inc1", "plane")
        inc1.set_axis_direction("-z")
        inc1.set_phase(0)
        inc1.set_by_min_max('x', 0, period)
        inc1.set_by_min_max('y', 0, period)
        inc1.set_by_value('z',25e-6)
        inc1.assemble()
        inc2 = fdtd.add_source("Inc2", "plane")
        inc2.set_axis_direction("-z")
        inc2.set_phase(-90)
        inc2.set_polarization(90)
        inc2.set_by_min_max('x', 0, period)
        inc2.set_by_min_max('y', 0, period)
        inc2.set_by_value('z',25e-6)
        inc2.assemble()

    ref = fdtd.add_monitor("ref", "dft")
    ref.set_by_min_max('x', 0, period)
    ref.set_by_min_max('y', 0, period)
    ref.set_by_value('z',28e-6)
    ref.keep_only_power()
    ref.assemble()

    fdtd.save()
    print(f"Running ...")
    fdtd.run()
    reflectance = fdtd.get_result("ref", "T")["T"]
    absorption = 1 - reflectance

    print(f"Closing ...")
    fdtd.close()

    return absorption
