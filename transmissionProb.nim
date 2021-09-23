
import math, strutils, algorithm, random, sequtils, os, strformat, tables
import json except `{}`

import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm

import ggplotnim



let siNfile03 = &"./resources/Si3N4Density=3.44Thickness=0.3microns" 
let siNfile015 = &"./resources/Si3N4Density=3.44Thickness=0.15microns" 
let siNfile01 = &"./resources/Si3N4Density=3.44Thickness=0.1microns"
let siFile = "./resources/SiDensity=2.33Thickness=200.microns"
let detectorFile = "./resources/transmission-argon-30mm-1050mbar-295K.dat"
let alFile002 = &"./resources/AlDensity=2.7Thickness=0.02microns"
let alFile0015 = &"./resources/AlDensity=2.7Thickness=0.015microns"

var dfTab = initTable[string, DataFrame]()
dfTab["siNfile0.3"] = readCsv(siNfile03, sep = ' ')
dfTab["siNfile0.15"] = readCsv(siNfile015, sep = ' ')
dfTab["siNfile0.1"] = readCsv(siNfile01, sep = ' ')
dfTab["siFile"] = readCsv(siFile, sep = ' ')
dfTab["detectorFile"] = readCsv(detectorFile, sep = ' ')
dfTab["alFile0.02"] = readCsv(alFile002, sep = ' ')
dfTab["alFile0.015"] = readCsv(alFile0015, sep = ' ')

let energies = dfTab["siNfile0.3"]["PhotonEnergy(eV)"].toTensor(float)
let transSiN03 = dfTab["siNfile0.3"]["Transmission"].toTensor(float)
let transSiN015 = dfTab["siNfile0.15"]["Transmission"].toTensor(float)
let transSiN01 = dfTab["siNfile0.1"]["Transmission"].toTensor(float)
let transSi = dfTab["siFile"]["Transmission"].toTensor(float)
let transDet = dfTab["detectorFile"]["Transmission"].toTensor(float)
let transAl002 = dfTab["alFile0.02"]["Transmission"].toTensor(float)
let transAl0015 = dfTab["alFile0.015"]["Transmission"].toTensor(float)

let dfTransProb = seqsToDf({"Photon energy [eV]": energies,
                              "Si3N4 0.3um": transSiN03,
                              "Si3N4 0.15um": transSiN015,
                              "Si3N4 0.1um": transSiN01,
                              "Si 200um": transSi,
                              #"TransAr": transDet,
                              "Al 0.02um": transAl002,
                              "Al 0.015um": transAl0015})#.mutate(f{"Ar 3cm" ~ 1.0 - `TransAr`}).arrange("Photon energy [eV]")
let dfTransProbNew = dfTransProb.gather(["Si3N4 0.3um", "Si3N4 0.15um", "Si3N4 0.1um", "Si 200um", #"Ar 3cm", 
                                         "Al 0.02um", "Al 0.015um"], key = "Material", value = "Transmission Probability")
echo dfTransProbNew
ggplot(dfTransProbNew) +
  geom_line(aes("Photon energy [eV]", "Transmission Probability", color = "Material")) +
  #ggtitle("The transmission probability for different window material with different thicknesses as well as the absorption probability of 3 cm Argon gas") +
  ggsave(&"out/TransProb.pdf")

let dfTransBest = seqsToDf({"Photon energy [eV]": energies,
                            "Si3N4": transSiN01,
                            "Al": transAl0015}).mutate(f{"SiNAl" ~ `Si3N4` * `Al`})

ggplot(dfTransBest) +
  geom_line(aes("Photon energy [eV]", "SiNAl")) +
  ggsave(&"out/TransProbBest.pdf")
