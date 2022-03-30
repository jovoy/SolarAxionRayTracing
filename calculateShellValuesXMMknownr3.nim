import math, strutils, algorithm, random, sequtils, os, strformat, tables
import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm

# Calculate the radii and angles of the XMM Telescope with the Information we have # all units are mm except angle

var
  ts: seq[float] #thickness of the shell from 0.47 mm to 1.07 mm innermost to outermost
  ds: seq[float] 
  r1: float
  r1j: float
  r1s: seq[float]
  r2s: seq[float]
  r4s: seq[float]
  r5s: seq[float]
  alpha: float
  alphas: seq[float]
  r2 : float
  r3: float
  r4: float
  r5: float
  xSep: float
  r1old: float
  focal = 7500.0
  lMirr = 300.0
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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #not useful #distance between the shell ranging from 1 mm to 5 mm from innermost to outermost shell
  allXsep = newSeq[float](R1.len)

for i in 1..58:  
  
  r3 = allR3[i-1]
  xSep = allXsep[i-1]
  alpha= radToDeg((arctan(r3/focal))) / 4.0
  r2 = r3 + 0.5 * xSep * tan(degToRad(alpha))
  r1 = r2 + lMirr * sin(degToRad(alpha))
  r4 = r3 - 0.5 * xSep * tan(3.0 * degToRad(alpha))
  r5 = r4 - lMirr * sin(3.0 * degToRad(alpha))
  alphas.add(round(alpha,3))
  r1s.add(round(r1,3))
  r2s.add(round(r2,3))
  r4s.add(round(r4,3))
  r5s.add(round(r5,3))

  

echo "R1 = ", r1s
echo "R2 = ", r2s
echo "R3 = ", allR3
echo "R4 = ", r4s
echo "R5 = ", r5s
echo "Thicknesses = ", allThickness
echo "f = ", focal
echo "xSep = ", allXsep
