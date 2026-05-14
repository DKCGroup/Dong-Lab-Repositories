import MOEMS
import numpy as np
from matplotlib import pyplot as plt
import os


os.listdir("./examples")
for _dir in os.listdir("./examples"):
    if os.path.isdir("./examples/" + _dir):
        print(f"Enter Path: ./examples/{_dir}")
        filepath = f"./examples/{_dir}/"
        filename = _dir
        print(f"Current Task: {_dir} Phase M LCP")
        LCP = MOEMS.MOEMS_generate_by_image(filepath, filename,
                                    source="LCP", phase="M",
                                    period=5e-6, h_Pt=50e-9,
                                    h_VO2=0.16e-6, h_Cr=30e-9,
                                    h_Au=90e-9, h_air=1e-6)
        print("Task Finished")
        print(f"Current Task: {_dir} Phase M RCP")
        RCP = MOEMS.MOEMS_generate_by_image(filepath, filename,
                                            source="RCP", phase="M",
                                            period=5e-6, h_Pt=50e-9,
                                            h_VO2=0.16e-6, h_Cr=30e-9,
                                            h_Au=90e-9, h_air=1e-6)
        print("Task Finished")
        plt.plot(np.linspace(10.1, 11.1, 200), LCP-RCP)
        plt.title("CD")
        plt.savefig(f"./examples/{_dir}/{_dir}_CD.png")
        plt.close()
        print("Exit Path: ./examples/{_dir}")
