import lumapi


MATERIAL = {
    "Ti": "Ti (Titanium) - Palik",
    "Ge": "Ge (Germanium) - Palik",
    "Al": "Al",
    "Cu": "Cu (Copper) - Palik",
    "Si": "Si (Silicon) - Palik",
    "W": "W (Tungsten) - Palik",
    "Au": "Au (Gold) - Palik",
    "Ag": "Ag (Silver) - Palik (1-10um)",
    "Pt": "Pt (Platinum) - Palik",
    "Cr": "Cr (Chromium) - Palik",
    "SiO2": "SiO2 (Glass) - Palik",
    "GaAs": "GaAs - Palik",
    "Si3N4": "Si3N4 (Silicon Nitride) - Kischkat",
    "Al2O3": "Al2O3 - Palik",
    "InAs": "InAs - Palik",
    "InP": "InP - Palik",
    "TiO2": "TiO2 (Titanium Dioxide) - Kischkat",
}


class Model:

    def __init__(self, filepath, filename, show=False):
        self.filepath = filepath
        self.filename = filename
        if show:
            self.model = lumapi.FDTD()
        else:
            self.model = lumapi.FDTD(hide=True)
        self.elements = {}
        self.groups = {}
        self.sims = {}
        self.sources = {}
        self.monitors = {}
        self.meshes = {}
        self.state = "layout"
        self.model.eval(f'cd("{filepath}");')
        self.model.save(f'{filename}.fsp')

    def switch_to_layout(self):
        if self.state == "layout":
            # print("You are already in layout mode!")
            pass
        else:
            self.model.eval('switchtolayout;\n')
            self.state = "layout"

    def import_material(self, filename):
        self.model.eval(f'importmaterialdb("{filename}");\n')

    def set_global_monitor(self, type="MIR"):
        if type == "MIR":
            code = 'setglobalmonitor("use wavelength spacing", true);\n'
            code += 'setglobalmonitor("use source limits", false);\n'
            code += 'setglobalmonitor("minimum wavelength", 8e-6);\n'
            code += 'setglobalmonitor("maximum wavelength", 13e-6);\n'
            code += 'setglobalmonitor("frequency points", 200);\n'
            self.model.eval(code)
        elif type == "10.6um":
            code = 'setglobalmonitor("use wavelength spacing", true);\n'
            code += 'setglobalmonitor("use source limits", false);\n'
            code += 'setglobalmonitor("minimum wavelength", 10.1e-6);\n'
            code += 'setglobalmonitor("maximum wavelength", 11.1e-6);\n'
            code += 'setglobalmonitor("frequency points", 200);\n'
            self.model.eval(code)
        else:
            print("Unknown Frequency Range.")

    def set_global_source(self, type="MIR"):
        if type == "MIR":
            code = 'setglobalsource("wavelength start", 8e-6);\n'
            code += 'setglobalsource("wavelength stop", 13e-6);\n'
            self.model.eval(code)
        elif type == "10.6um":
            code = 'setglobalsource("wavelength start", 10.1e-6);\n'
            code += 'setglobalsource("wavelength stop", 11.1e-6);\n'
            self.model.eval(code)
        else:
            print("Unknown Frequency Range.")

    def add_structure_group(self, name):
        group = StructureGroup(self.model, name)
        self.groups[name] = group
        return group

    def add_rect(self, name):
        rect = Element(self.model, name, f'addrect;\nset("name", "{name}");\n')
        self.elements[name] = rect
        return rect

    def import_gds(self, name, filename, cellname, layer, material, zmin, zmax):
        self.model.eval(f'gdsimport("{filename}.gds", "{cellname}", "{layer}", "{MATERIAL.get(material, material)}", {zmin}, {zmax});')
        group = StructureGroup(self.model, name, exists=True)
        self.groups[name] = group
        return group

    def add_fdtd(self, name="fdtd"):
        sim = Sim(self.model, name)
        self.sims[name] = sim
        return sim

    def add_mesh(self, name):
        mesh = Mesh(self.model, name)
        self.meshes[name] = mesh
        return mesh

    def add_source(self, name, type):
        source = Source(self.model, name, type)
        self.sources[name] = source
        return source

    def add_monitor(self, name, type):
        monitor = Monitor(self.model, name, type)
        self.monitors[name] = monitor
        return monitor

    def save(self):
        self.model.save(f'{self.filename}.fsp')

    def run(self):
        self.model.save(f'{self.filename}.fsp')
        self.model.run()
        self.state = "analysis"

    def close(self):
        self.model.close()

    def get_result(self, monitor, var):
        return self.model.getresult(monitor, var)


class StructureGroup:
    def __init__(self, model, name, exists=False):
        self.model = model
        self.name = name
        self.elements = {}
        if not exists:
            self.code = f'addstructuregroup;\nset("name", "{name}");\nset("x", 0);\nset("y", 0);\n'
        else:
            self.code = ''

    def assemble(self):
        self.model.eval(self.code)

    def add_rect(self, name):
        rect = Element(self.model, name, f'addrect;\naddtogroup("{self.name}");\nset("name", "{name}");\n')
        self.elements[name] = rect
        return rect


class Element:

    def __init__(self, model, name, code):
        self.model = model
        self.name = name
        self.code = code

    def assemble(self):
        self.model.eval(self.code)

    def set_by_min_max(self, var, min, max):
        self.code += f'set("{var} min", {min});\nset("{var} max", {max});\n'

    def set_by_center_span(self, var, center, span):
        self.code = self.code  # TBD

    def set_material(self, mat):
        if mat in MATERIAL.keys():
            self.code += f'set("material", "{MATERIAL[mat]}");\n'
        else:
            self.code += f'set("material", "{mat}");\n'

    def set_mesh_order(self, order):
        self.code += f'set("override mesh order from material database", true);\n'
        self.code += f'set("mesh order", {order});\n'


class Sim:

    def __init__(self, model, name):
        self.model = model
        self.name = name
        self.code = 'addfdtd;\n'

    def assemble(self):
        self.model.eval(self.code)

    def set_by_min_max(self, var, min, max):
        self.code += f'set("{var} min", {min});\nset("{var} max", {max});\n'

    def set_boundary_condition(self, var, type):
        self.code += f'set("{var} min bc", "{type}");\n'
        self.code += f'set("{var} max bc", "{type}");\n'

    def set_simulation_time(self, time=5000):
        self.code += f'set("simulation time", {time}*1e-15);\n'

    def set_auto_shutoff(self, shutoff=1e-5):
        self.code += f'set("auto shutoff min", {shutoff});\n'

    def set_mesh_type(self, type="custom"):
        if type == "custom":
            self.code += f'set("mesh type", "custom non-uniform");\n'
        elif type == "uniform":
            self.code += f'set("mesh type", "uniform");\n'
        else:
            self.code += f'set("mesh type", "auto non-uniform");\n'

    def set_mesh_definition(self, axis, type="mesh step"):
        if type == "wavelength":
            self.code += f'set("define {axis} mesh by", "mesh cells per wavelength");\n'
        if type == "mesh step":
            self.code += f'set("define {axis} mesh by", "maximum mesh step");\n'
        if type == "hybrid":
            self.code += f'set("define {axis} mesh by", "max mesh step and mesh cells per wavelength");\n'
        if type == "number":
            self.code += f'set("define {axis} mesh by", "number of mesh cells");\n'

    def set_maximum_mesh_step(self, axis, step=0.1):
        self.code += f'set("d{axis}", {step}*1e-6);\n'

    def set_mesh_cell_per_wavelength(self, cells=10):
        self.code += f'set("mesh cells per wavelength", {cells});\n'

    def set_dt_stability_factor(self, factor=0.95):
        self.code += f'set("dt stability factor", {factor});\n'

    def set_minimum_mesh_step(self, step=0.001):
        self.code += f'set("min mesh step", {step}*1e-6);\n'


class Source:

    def __init__(self, model, name, type):
        self.model = model
        self.name = name
        if type == 'plane':
            self.code = f'addplane;\n'
            # TBD
        self.code += f'set("name", "{name}");\n'

    def assemble(self):
        self.model.eval(self.code)

    def set_axis_direction(self, axis):
        if axis == "-z":
            self.code += f'set("injection axis","z");\n'
            self.code += f'set("direction","backward");\n'
            # TBD

    def set_phase(self, phase):
        self.code += f'set("phase", {phase});\n'

    def set_polarization(self, pol):
        self.code += f'set("polarization angle", {pol});\n'
        if pol == "TM":
            self.code += f'set("polarization angle", 90);\n'

    def set_by_min_max(self, var, min, max):
        self.code += f'set("{var} min", {min});\nset("{var} max", {max});\n'

    def set_by_value(self, var, value):
        self.code += f'set("{var}", {value});\n'

    def set_wavelength(self, min, max):
        self.code += f'set("wavelength start", {min});\n'
        self.code += f'set("wavelength stop", {max});\n'


class Monitor:

    def __init__(self, model, name, type):
        self.model = model
        self.name = name
        if type == 'dft':
            self.code = f'adddftmonitor;\n'
            # TBD
        self.code += f'set("name", "{name}");\n'

    def assemble(self):
        self.model.eval(self.code)

    def set_axis_direction(self, axis):
        if axis == "-z":
            self.code += 'set("monitor type",7);\n'
        if axis == "-x":
            self.code += 'set("monitor type",5);\n'

    def set_by_min_max(self, var, min, max):
        self.code += f'set("{var} min", {min});\nset("{var} max", {max});\n'

    def set_by_value(self, var, value):
        self.code += f'set("{var}", {value});\n'

    def set_freq_points(self, num):
        self.code += f'setglobalmonitor("frequency points",{num});'

    def keep_only_power(self):
        self.code += f'set("output Ex", false);\n'
        self.code += f'set("output Ey", false);\n'
        self.code += f'set("output Ez", false);\n'
        self.code += f'set("output Hx", false);\n'
        self.code += f'set("output Hy", false);\n'
        self.code += f'set("output Hz", false);\n'


class Mesh:
    def __init__(self, model, name):
        self.model = model
        self.name = name
        self.code = f'addmesh;\n'
        self.code += f'set("name", "{name}");\n'

    def assemble(self):
        self.model.eval(self.code)

    def set_by_min_max(self, var, min, max):
        self.code += f'set("{var} min", {min});\nset("{var} max", {max});\n'

    def set_by_value(self, var, value):
        self.code += f'set("{var}", {value});\n'

    def set_override_direction(self, dir):
        if dir == "x":
            self.code += f'set("override x mesh", 1);\n'
            self.code += f'set("override y mesh", 0);\n'
            self.code += f'set("override z mesh", 0);\n'
        if dir == "y":
            self.code += f'set("override y mesh", 1);\n'
            self.code += f'set("override x mesh", 0);\n'
            self.code += f'set("override z mesh", 0);\n'
        if dir == "z":
            self.code += f'set("override z mesh", 1);\n'
            self.code += f'set("override x mesh", 0);\n'
            self.code += f'set("override y mesh", 0);\n'
        if dir == "xy":
            self.code += f'set("override x mesh", 1);\n'
            self.code += f'set("override y mesh", 1);\n'
            self.code += f'set("override z mesh", 0);\n'
        if dir == "yz":
            self.code += f'set("override y mesh", 1);\n'
            self.code += f'set("override z mesh", 1);\n'
            self.code += f'set("override x mesh", 0);\n'
        if dir == "xz":
            self.code += f'set("override x mesh", 1);\n'
            self.code += f'set("override z mesh", 1);\n'
            self.code += f'set("override y mesh", 0);\n'
        if dir == "xyz":
            self.code += f'set("override x mesh", 1);\n'
            self.code += f'set("override y mesh", 1);\n'
            self.code += f'set("override z mesh", 1);\n'
    
    def set_override_mesh_step(self, axis, step=0.1):
        self.code += f'set("set maximum mesh step",1);\n'
        self.code += f'set("d{axis}", {step}*1e-6);\n'

