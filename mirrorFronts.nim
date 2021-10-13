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
  allXsep = @[4.0, 4.171, 4.140, 4.221, 4.190, 4.228, 4.245, 4.288, 4.284,
          4.306, 4.324, 4.373, 4.387, 4.403, 4.481]
  rBore = 21.5
  dGlas = 0.2
  d = allR1[0] + rBore 
  xGR = 2.0
  totalDisk = PI * rBore * rBore

var totalArea: float

for i in 1..14:
  let 
    rBot = allR1[i-1] + dGlas
    aI = rBore * rBore * arccos((d * d + rBore * rBore - allR1[i] * allR1[i]) / (2.0 * d * rBore)) +
        allR1[i] * allR1[i] * arccos((d * d - rBore * rBore + allR1[i] * allR1[i]) / (2.0 * d * allR1[i] * allR1[i])) -
        0.5 * sqrt((- d + rBore + allR1[i]) * (d + rBore - allR1[i])) *
        sqrt((d - rBore + allR1[i]) * (d + rBore + allR1[i])) 
    aI_1 = rBore * rBore * arccos((d * d + rBore * rBore - rBot * rBot) / (2.0 * d * rBore)) +
        rBot * rBot * arccos((d * d - rBore * rBore + rBot * rBot) / (2.0 * d * rBot * rBot)) -
        0.5 * sqrt((- d + rBore + rBot) * (d + rBore - rBot)) *
        sqrt((d - rBore + rBot) * (d + rBore + rBot))
    aCSI = aI - aI_1 - (allR1[i] - allR1[i-1]) * xGR
  echo aI
  echo aI_1
  echo aCSI
  totalArea += aCSI

echo totalArea
echo totalDisk