import nimpy, ggplotnim, nimhdf5, unchained
import scinim / numpyarrays
from std / os import `/`

#[
Simple tool to generate an HDF5 file for the reflectivities of all the
different layers of the LLNL telescope.

It requires the Python library `darpanx`.

Run the script without arguments to generate the file that we need for us.

To generate the file as required by REST, run instead:

```sh
./llnl_layer_reflectivity \
    --numEnergy 500 \
    --angleMax 9.0.° \
    --numAngle 901 \
    --outfile llnl_layer_reflectivities_rest.h5 \
    --outpath /path/to/foo
```
]#

let drp = pyImport("darpanx")

let m1 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 115.0,
                        D_max = 225.0,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=2,
                        SigmaValues=[1.0])

let m2 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 70.0,
                        D_max = 190.0,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=3,
                        SigmaValues=[1.0])

let m3 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 55.0,
                        D_max = 160.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=4,
                        SigmaValues=[1.0])

let m4 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 50.0,
                        D_max = 140.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=5,
                        SigmaValues=[1.0])

let ms = [m1, m2, m3, m4]

proc main(energyMin = 0.03.keV,
          energyMax = 15.0.keV,
          numEnergy = 1000,
          angleMin = 0.0.°,
          angleMax = 1.5.°,
          numAngle = 1000,
          outfile = "llnl_layer_reflectivities.h5",
          outpath = "../resources/") =
  let Energy = linspace(energyMin.float, energyMax.float, numEnergy) # scan in a range wider than we use in the raytracer
  let θs = linspace(angleMin.float, angleMax.float, numAngle) # [0.5] # the angle under which we check

  var h5f = H5open(outpath / outfile, "rw")
  var eDset = h5f.create_dataset("Energy", (Energy.len, 1), dtype = float64)
  var aDset = h5f.create_dataset("Angles", (θs.len, 1), dtype = float64)
  eDset[eDset.all] = Energy
  aDset[aDset.all] = θs

  let Enp = Energy.toTensor.toNdArray
  for i, m in ms:
    var df = newDataFrame()
    var refl = zeros[float]([θs.len, Energy.len])
    var ϑidx = 0
    for ϑ in θs:
      # compute the reflectivity
      discard m.get_optical_func(Theta = [ϑ], Energy = Enp)
      let reflect = toTensor[float64](m.Ra)
      refl[ϑidx, _] = reflect.unsqueeze(axis = 0)
      let dfLoc = seqsToDf({"Energy" : Energy, "Theta" : ϑ, "Ref" : reflect})
      df.add dfLoc
      inc ϑidx
    ggplot(df, aes("Energy", "Theta", fill = "Ref")) +
      geom_raster() +
      ggtitle("Reflectivity of LLNL layer type " & $i) +
      ggsave("/tmp/layer" & $i & ".pdf")
    var rDset = h5f.create_dataset("Reflectivity" & $i, (θs.len, Energy.len),
                                   dtype = float64)
    rDset.unsafeWrite(refl.dataArray, refl.size)

when isMainModule:
  import cligen
  import unchained / cligenParseUnits
  dispatch main
