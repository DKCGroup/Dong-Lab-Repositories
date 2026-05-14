import TARC
import numpy as np
from matplotlib import pyplot as plt
import os

os.listdir("./examples")
for _dir in os.listdir("./examples"):
    if os.path.isdir("./examples/" + _dir):
        print(f"Enter Path: ./examples/{_dir}")
        filepath = f"./examples/{_dir}/"
        filename = _dir
        print(f"Current Task: {_dir} Phase M")
        M = TARC.TARC_generate_by_image(filepath, filename,phase="M",
                                    period=5e-6, height=2.5e-6)
        print("Task Finished")
        print(f"Current Task: {_dir} Phase I")
        I = TARC.TARC_generate_by_image(filepath, filename,phase="I",
                                    period=5e-6, height=2.5e-6)
        print("Task Finished")
        plt.plot(np.linspace(8, 13, 200), I, label="I")
        plt.plot(np.linspace(8, 13, 200), M, label="M")
        plt.ylim(0, 1)
        plt.ylabel("Absorption")
        plt.legend()
        plt.savefig(f"./examples/{_dir}/{_dir}_spectrum.png")
        plt.close()
        print("Exit Path: ./examples/{_dir}")
