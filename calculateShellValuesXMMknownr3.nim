import math, strutils, algorithm, random, sequtils, os, strformat, tables
import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm

# Calculate the radii and angles of the XMM Telescope with the Information we have # all units are mm except angle

var
  t: float
  ts: seq[float] #thickness of the shell from 0.47 mm to 1.07 mm innermost to outermost
  d: float
  ds: seq[float] #distance between the shell ranging from 1 mm to 5 mm from innermost to outermost shell
  r1 = 153.0
  r1j: float
  r1s: seq[float]
  alpha: float
  alphas: seq[float]
  x = 0
  xsep: seq[float]
  r3: float
  r5: float
  r1old: float
  focal = 7500.0
let 
  allR3 = @[151.61, 153.88, 156.17, 158.48, 160.82, 163.18, 165.57, 167.98, 170.42, 172.88, 175.37, 177.88, 180.42, 183.14, 185.89, 188.67, 191.48, 
      194.32, 197.19, 200.09, 203.02, 206.03, 209.07, 212.14, 215.24, 218.37, 221.54, 224.74, 227.97, 231.24, 234.54, 237.87, 241.24, 244.85, 
      248.5, 252.19, 255.92, 259.68, 263.48, 267.32, 271.2, 275.12, 279.08, 283.09, 287.14, 291.38, 295.72, 300.11, 304.54, 309.02, 313.54, 
      318.11, 322.73, 327.4, 332.12, 336.88, 341.69, 346.55] #these are the real values but R3
  allThickness = @[0.468, 0.475, 0.482, 0.490, 0.497, 0.504, 0.511, 0.519, 0.526, 0.534, 0.542, 0.549, 0.557, 0.566, 0.574, 0.583, 0.591, 0.600, 0.609, 0.618, 
      0.627, 0.636, 0.646, 0.655, 0.665, 0.675, 0.684, 0.694, 0.704, 0.714, 0.724, 0.735, 0.745, 0.756, 0.768, 0.779, 0.790, 0.802, 0.814, 0.826, 0.838, 0.850, 
      0.862, 0.874, 0.887, 0.900, 0.913, 0.927, 0.941, 0.955, 0.968, 0.983, 0.997, 1.011, 1.026, 1.041, 1.055, 1.070] #these are only theoretical values but official ones
  allDistances = @[0.294, 0.284, 0.274, 0.273, 0.263, 0.263, 0.252, 0.250, 0.239, 0.236, 0.223, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.465, 0.458, 0.450, 0.442, 0.424, 0.415, 0.405, 0.395, 0.385, 0.374, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

for i in 1..58:  
  ts = allThickness  
  ds = allDistances
  
  
  r3 = allR3[i-1]
  alpha= radToDeg((arctan(r3/focal))) / 4.0
  r1 = r3 + 300.0 * sin(degToRad(alpha))
  r5 = r3 - 300.0 * sin(3.0*degToRad(alpha))
  echo "r5 ", r5
  alphas.add(round(alpha,3))
  xsep.add(0.0)
  if ds[i-1] == 0.0:
    ds[i-1] = 0.5
  echo ds[i-1]
  r1s.add(round(r1,3))
  r1old = r1 + ts[i-1] + ds[i-1]
  

echo x
echo r1s
echo xsep
echo alphas
