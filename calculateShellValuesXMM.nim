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

for i in 1..58:
  t = 0.47+(i.float-1.0)*((1.07-0.47)/57.0)
  ts.add(t)
  d = 1.0+(i.float-1.0)*((5.0-1.0)/57.0)
  ds.add(d)
  
  r1s.add(round(r1,3))
  
  for j in 1..100:
    r1j = 300.0*sin(degToRad(j.float*0.01))+7500*tan(4.0*degToRad(j.float*0.01))
    if r1j < r1+2.7 and r1j > r1-2.7:
      echo r1
      echo r1j
      echo j.float*0.01
      alpha = j.float*0.01
      x+=1

  alphas.add(round(alpha,3))
  xsep.add(0.0)
  r1 = r1 + ts[i-1] + ds[i-1]
  

echo x
echo r1s
echo xsep
echo alphas
