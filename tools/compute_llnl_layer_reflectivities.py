# Needs the DarpanX python library
# https://github.com/biswajitmb/DarpanX

# Based on the DarpanX example 3

# Computes the reflectivities of the different kinds of layers of the
# LLNL telescope. See the DTU PhD thesis to understand where these
# parameters here come from

import darpanx as drp
import numpy as np

m1 = drp.Multilayer(MultilayerType="DepthGraded",
                    SubstrateMaterial="SiO2",
                    D_min = 11.5,
                    D_max = 22.5,
                    Gamma = 0.45,
                    C = 1.0,
                    LayerMaterial=["C", "Pt"],
                    Repetition=2,
                    SigmaValues=[1.0])

m2 = drp.Multilayer(MultilayerType="DepthGraded",
                    SubstrateMaterial="SiO2",
                    D_min = 7.0,
                    D_max = 19.0,
                    Gamma = 0.45,
                    C = 1.0,
                    LayerMaterial=["C", "Pt"],
                    Repetition=3,
                    SigmaValues=[1.0])

m3 = drp.Multilayer(MultilayerType="DepthGraded",
                    SubstrateMaterial="SiO2",
                    D_min = 5.5,
                    D_max = 16.0,
                    Gamma = 0.4,
                    C = 1.0,
                    LayerMaterial=["C", "Pt"],
                    Repetition=4,
                    SigmaValues=[1.0])

m4 = drp.Multilayer(MultilayerType="DepthGraded",
                    SubstrateMaterial="SiO2",
                    D_min = 5.0,
                    D_max = 14.0,
                    Gamma = 0.4,
                    C = 1.0,
                    LayerMaterial=["C", "Pt"],
                    Repetition=5,
                    SigmaValues=[1.0])


ms = [m1, m2, m3, m4]

Energy = np.linspace(0.01, 15.0, 1000) # scan in a range wider than we use in the raytracer
Theta = [0.5] # the angle under which we check

from pathlib import Path
Path("plots").mkdir(exist_ok = True)
for i, m in enumerate(ms):
    m.get_optical_func(Theta = Theta, Energy = Energy)
    m.plot(ylog = "no", Comp = ["Ra"], OutFile = "plots/Pt_SiC_" + str(i), Struc = 'yes')
