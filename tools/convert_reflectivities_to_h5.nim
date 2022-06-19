import std / [strformat, strutils, os, sugar, strscans, algorithm]
import nimhdf5, arraymancer, ggplotnim

## A simple tool to convert the gold reflectivity files downloaded from henke.gov
## into a single H5 file for convenience.

# Constants for the output of the downloaded files
const outpath = "../resources/henke_download/"

# collect it
let goldFiles = collect(newSeq):
  for f in walkFiles(outpath & "*degGold0.25microns.csv"):
    let (success, angle) = scanTuple(f.dup(removePrefix(outpath)), "$fdeg")
    if not success: raise newException(IOError, "Could not parse input gold file: " & $f)
    (f, angle)

# sort it
let sortedFiles = goldFiles.sortedByIt(it[1])

# read it
var θs = newSeq[float]()
var goldReflect = newTensor[float](0)
let Energy = linspace(0.03, 15.0, 1000) # scan in a range wider than we use in the raytracer
var df = newDataFrame()
for i, (f, angle) in sortedFiles:
  θs.add angle
  var dfLoc = readCsv(f, sep = ' ', header = "#")
  if goldReflect.size == 0: # means first iteration, don't know # elements in gold files
    goldReflect = newTensor[float]([goldFiles.len, dfLoc.len])
  let reflect = dfLoc["Reflectivity", float]
  goldReflect[i, _] = reflect.unsqueeze(0)
  df.add seqsToDf({"Energy" : Energy, "Theta" : angle, "Ref" : reflect})

# plot it
ggplot(df, aes("Energy", "Theta", fill = "Ref")) +
  geom_raster() +
  ggtitle("Reflectivity of Gold 250nm") +
  ggsave("/tmp/gold_0.25microns_reflectivity.pdf")
# write it
var h5f = H5open("../resources/gold_0.25microns_reflectivities.h5", "rw")
var eDset = h5f.create_dataset("Energy", (Energy.len, 1), dtype = float64)
var aDset = h5f.create_dataset("Angles", (θs.len, 1), dtype = float64)
eDset[eDset.all] = Energy
aDset[aDset.all] = θs

var rDset = h5f.create_dataset("Reflectivity", (Energy.len, θs.len),
                               dtype = float64)
rDset.unsafeWrite(goldReflect.dataArray, goldReflect.size)
