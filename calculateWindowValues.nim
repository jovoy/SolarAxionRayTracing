import math, strutils, algorithm, random, sequtils, os, strformat, tables
import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim
import glm

let 
  radiusWindow = 7.0 #4.0
  openAperatureRatio = 0.838 #0.95
  totalArea = radiusWindow * radiusWindow * PI
  areaOfStrips = totalArea * (1.0 - openAperatureRatio)
  numberOfStrips = 4 #20 #has to be even so there is none in the middle
  dAndwPerStrip = radiusWindow * 2.0 / (numberOfStrips.float + 1.0)  #width and distance between strips per strip; 
                                                                     #the width on both sides is a whole width 
                                                                     #(Don't know what to do about the additional string width) 
                                                                     #not important at high strip number but at low

var
  lengthStrip: float
  lengthAllStrips: float

for i in 0..(numberOfStrips/2).round.int - 1:
  lengthStrip = sqrt(radiusWindow * radiusWindow - (i.float * dAndwPerStrip + 0.5 * dAndwPerStrip) * (i.float * dAndwPerStrip + 0.5 * dAndwPerStrip)) * 2.0
  lengthAllStrips += lengthStrip
  echo lengthStrip

lengthAllStrips = lengthAllStrips * 2.0

let 
  widthStrips = areaOfStrips / lengthAllStrips
  distStrips = dAndwPerStrip - widthStrips
echo widthStrips
echo distStrips
