import nimpy, ggplotnim, nimhdf5
import scinim / numpyarrays

#let np = pyImport("numpy")
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

let Energy = linspace(0.01, 15.0, 1000) # scan in a range wider than we use in the raytracer
let ϑs = linspace(0.0, 1.5, 1000) # [0.5] # the angle under which we check

var h5f = H5open("/tmp/data.h5", "rw")
var eDset = h5f.create_dataset("Energy", (Energy.len, 1), dtype = float64)
var aDset = h5f.create_dataset("Angles", (ϑs.len, 1), dtype = float64)
eDset[eDset.all] = Energy
aDset[aDset.all] = ϑs

let Enp = Energy.toTensor.toNdArray
for i, m in ms:
  var df = newDataFrame()
  var refl = zeros[float]([Energy.len, ϑs.len])
  var ϑidx = 0
  for ϑ in ϑs:
    # compute the reflectivity
    discard m.get_optical_func(Theta = [ϑ], Energy = Enp)
    let reflect = toTensor[float64](m.Ra)
    refl[ϑidx, _] = reflect
    let dfLoc = seqsToDf({"Energy" : Energy, "Theta" : ϑ, "Ref" : reflect})
    df.add dfLoc
    inc ϑidx
  echo df
  #df.writeCsv("/tmp/reflectivity_angles.csv")
  ggplot(df, aes("Energy", "Theta", fill = "Ref")) +
    geom_raster() +
    ggtitle("Reflectivity of LLNL layer type " & $i) +
    ggsave("/tmp/layer" & $i & ".pdf")
  var rDset = h5f.create_dataset("Reflectivity" & $i, (Energy.len, ϑs.len), dtype = float64)
  rDset.unsafeWrite(refl.dataArray, refl.size)
