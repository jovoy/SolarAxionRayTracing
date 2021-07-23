# stdlib
import math, strutils, algorithm, random, sequtils, os, strformat, tables
import json except `{}`

# TODO: axionMass will be minified and the important stuff extracted
# import axionMass/axionMass

# nimble
import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm
# import plotly
import ggplotnim


let 
  allR3 = @[151.61, 153.88, 156.17, 158.48, 160.82, 163.18, 165.57, 167.98, 170.42, 172.88, 175.37, 177.88, 180.42, 183.14, 185.89, 188.67, 191.48, 
      194.32, 197.19, 200.09, 203.02, 206.03, 209.07, 212.14, 215.24, 218.37, 221.54, 224.74, 227.97, 231.24, 234.54, 237.87, 241.24, 244.85, 
      248.5, 252.19, 255.92, 259.68, 263.48, 267.32, 271.2, 275.12, 279.08, 283.09, 287.14, 291.38, 295.72, 300.11, 304.54, 309.02, 313.54, 
      318.11, 322.73, 327.4, 332.12, 336.88, 341.69, 346.55] #these are the real values but R3
  allR1 = @[153.126, 155.419, 157.731, 160.065, 162.428, 164.812, 167.225, 169.66, 172.124, 174.608, 177.123, 179.658, 182.224, 184.971, 187.749, 190.556, 
      193.394, 196.263, 199.161, 202.09, 205.05, 208.09, 211.16, 214.261, 217.392, 220.553, 223.755, 226.987, 230.249, 233.552, 236.885, 240.248, 243.652, 
      247.298, 250.984, 254.711, 258.478, 262.276, 266.114, 269.992, 273.911, 277.87, 281.869, 285.92, 290.01, 294.292, 298.676, 303.109, 307.584, 312.108, 
      316.674, 321.289, 325.955, 330.672, 335.439, 340.246, 345.104, 350.013] ## the radii of the shells closest to the magnet, now correct
  allThickness = @[0.468, 0.475, 0.482, 0.490, 0.497, 0.504, 0.511, 0.519, 0.526, 0.534, 0.542, 0.549, 0.557, 0.566, 0.574, 0.583, 0.591, 0.600, 0.609, 0.618, 
      0.627, 0.636, 0.646, 0.655, 0.665, 0.675, 0.684, 0.694, 0.704, 0.714, 0.724, 0.735, 0.745, 0.756, 0.768, 0.779, 0.790, 0.802, 0.814, 0.826, 0.838, 0.850, 
      0.862, 0.874, 0.887, 0.900, 0.913, 0.927, 0.941, 0.955, 0.968, 0.983, 0.997, 1.011, 1.026, 1.041, 1.055, 1.070]
  coatingMaterial = "Gold"
  coatingDensity = 19.302 #g/cm³
  coatingThickness = 0.00025 #mm = 250 nm
  totalArea = 3848.3375 #cm^2
  
var 
  areaMirrorFronts:float
  goldfile: string
  alpha: float
  alpha1: float
  dfTable = initTable[string, DataFrame]()
  energies: seq[float]
  reflectivity: seq[float]
  reflectionProbE: float
  reflectionProb: float
  reflectionEnergy: float

for i in 17..80:
  alpha = (i.float * 0.01).round(2)
  goldfile = fmt"./resources/reflectivity/{alpha:4.2f}degGold0.25microns"
  dfTable[fmt"goldfile{alpha:4.2f}"] = readCsv(goldfile, sep = ' ')


for l in 0..499:
  reflectionProbE = 0.0
  for j in 17..80:    
    alpha1 = (j.float * 0.01).round(2)
      #alpha2 = angle2[1].round(2)
    
    reflectionProb = dfTable[fmt"goldfile{alpha1:4.2f}"]["Reflectivity"].toTensor(float)[l]
    reflectionEnergy = dfTable[fmt"goldfile{alpha1:4.2f}"]["PhotonEnergy(eV)"].toTensor(float)[l]
    reflectionProbE += reflectionProb * 0.01
  reflectivity.add(reflectionProbE)
  energies.add(reflectionEnergy)

for k in 0..57:
  var areaRing = ((allR1[k] + allThickness[k]) * (allR1[k] + allThickness[k]) * PI) - (allR1[k] * allR1[k] * PI)
  areaMirrorFronts += areaRing * 0.01 #in cm^2
echo areaMirrorFronts

######## This part is about the spider and all of the values are estimated ##################


let dfRef = seqsToDf({"Reflectivity": reflectivity,
                      "Photon Energy[eV]": energies})
ggplot(dfRef.mutate(f{"EffArea" ~ `Reflectivity` * (totalArea - areaMirrorFronts)}),
        aes("Photon Energy[eV]", "EffArea")) +
  geom_line() +
  #xlim(100.0, 10000.0) +
  scale_x_log10() +
  #scale_y_log10() +
  ggtitle("The reflectivity od the XMM telescope depending on the photon energy") +
  ggsave(&"out/Reflectivity.pdf")
