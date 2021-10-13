import math, strutils, algorithm, random, sequtils, os, strformat, tables
import json except `{}`

import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm

import ggplotnim

let
  allR1 = @[60.7095, 63.006, 65.606, 68.305, 71.105, 74.011, 77.027, 80.157,
          83.405, 86.775, 90.272, 93.902, 97.668, 101.576, 105.632]
  rBore = 21.5
  dGlas = 0.2
  d = allR1[0] + rBore 
  xGR = 2.0
  totalDisk = PI * rBore * rBore

var 
  totalLostArea: float
  totalMirrorFronts: float
  totalGraphiteBlocks: float

for i in 1..<14:
  let 
    rBot = allR1[i] + dGlas
    aI = rBore * rBore * arccos((d * d + rBore * rBore - allR1[i] * allR1[i]) / (2.0 * d * rBore)) +
        allR1[i] * allR1[i] * arccos((d * d - rBore * rBore + allR1[i] * allR1[i]) / (2.0 * d * allR1[i])) -
        0.5 * sqrt((- d + rBore + allR1[i]) * (d + rBore - allR1[i])) *
        sqrt((d - rBore + allR1[i]) * (d + rBore + allR1[i])) 
    aI_gl = rBore * rBore * arccos((d * d + rBore * rBore - rBot * rBot) / (2.0 * d * rBore)) +
        rBot * rBot * arccos((d * d - rBore * rBore + rBot * rBot) / (2.0 * d * rBot)) -
        0.5 * sqrt((- d + rBore + rBot) * (d + rBore - rBot)) *
        sqrt((d - rBore + rBot) * (d + rBore + rBot))
    aCSI = aI_gl - aI + (allR1[i] - (allR1[i-1] + dGlas)) * xGR
  echo aI
  echo aI_gl
  echo aI_gl - aI
  totalMirrorFronts += aI_gl - aI
  totalGraphiteBlocks += (allR1[i] - (allR1[i-1] + dGlas)) * xGR
  
totalGraphiteBlocks += (d + rBore - (allR1[13] + dGlas)) * xGR #seems correct
totalLostArea = totalMirrorFronts + totalGraphiteBlocks
echo "results"
echo totalMirrorFronts
echo totalGraphiteBlocks
echo totalLostArea 
echo totalDisk