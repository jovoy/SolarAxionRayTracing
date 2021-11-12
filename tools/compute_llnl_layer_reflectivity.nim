## Needs the DarpanX python library
## https://github.com/biswajitmb/DarpanX
#
## Based on the DarpanX example 3
#
## Computes the reflectivities of the different kinds of layers of the
## LLNL telescope. See the DTU PhD thesis to understand where these
## parameters here come from
##
## Nim version using nimpy

import std / os
import nimpy, ggplotnim

let np = pyImport("numpy")
let drp = pyImport("darpanx")

let m1 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 11.5,
                        D_max = 22.5,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["Pt", "C"],
                        Repetition=2,
                        SigmaValues=[1.0])

let m2 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 7.0,
                        D_max = 19.0,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["Pt", "C"],
                        Repetition=3,
                        SigmaValues=[1.0])

let m3 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 5.5,
                        D_max = 16.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["Pt", "C"],
                        Repetition=4,
                        SigmaValues=[1.0])

let m4 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 5.0,
                        D_max = 14.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["Pt", "C"],
                        Repetition=5,
                        SigmaValues=[1.0])


let ms = [m1, m2, m3, m4]

let Energy = np.linspace(0.01, 15.0, 1000) # scan in a range wider than we use in the raytracer
let ϑ = [0.5] # the angle under which we check

createDir("plots")
for i, m in ms:
  discard m.get_optical_func(Theta = ϑ, Energy = Energy)
  discard m.plot(ylog = "no", Comp = ["Ra"], OutFile = "plots/Pt_SiC_" & $i, Struc = "yes")
