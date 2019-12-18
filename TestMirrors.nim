
import helpers/utils
import unittest

import glm/vec
import math
import strutils
import algorithm
import sequtils, os, strutils
include raytracer2018



#[proc findPosXRT(pointXRT : Vec3, pointCB : Vec3, r1 : float, r2 : float, angle : float, lMirror : float, distMirr : float, uncer : float, sMin : float, sMax : float, pointOrAngle : string) : Vec3 =

  ##this is to find the position the ray hits the mirror shell of r1. it is after transforming the ray into a cylindrical coordinate system, that has the middle and the beginning of the mirror "cylinders" as its origin and is in cylindrical coordinates
  
  var
    point = pointCB
    s : float
    term : float
    sMinHigh = (sMin * 1000000.0).int
    sMaxHigh = (sMax * 1000000.0).int
    pointMirror = vec3(0.0)
  let direc = pointXRT - pointCB
  
    
  for i in sMinHigh..sMaxHigh:
    s = i.float / 1000000.0
    term = sqrt((point[0] + s * direc[0]) * (point[0] + s * direc[0]) + (point[1] + s * direc[1]) *  (point[1] + s * direc[1])) - ((r2 - r1) * (point[2] + s * direc[2] - distMirr) / (cos(angle) * lMirror))
      
    if r1 < term + uncer and r1 > term - uncer:
      pointMirror = point + s * direc ## sometimes there are 3 different s for which this works, in that case the one with the highest value is taken
  echo "now"
  echo pointMirror
  var normalVec = vec3(0.0)
  normalVec[0] = pointMirror[0]
  normalVec[1] = pointMirror[1]
  normalVec[2] = tan(angle) * sqrt(pointMirror[0] * pointMirror[0] + pointMirror[1] * pointMirror[1])
  var vectorBeforeMirror =  pointXRT - pointCB
  var vectorAxis = vec3(0.0) ## this is the vector product of the normal vector on pointMirror pointing in the direction of the radius of the cylinder and the vector of the ray
  vectorAxis[0] = normalVec[1] * vectorBeforeMirror[2] - normalVec[2] * vectorBeforeMirror[1] #+ pointMirror1ZylKart[0]
  vectorAxis[1] = normalVec[2] * vectorBeforeMirror[0] - normalVec[0] * vectorBeforeMirror[2] #+ pointMirror1ZylKart[1]
  vectorAxis[2] = normalVec[0] * vectorBeforeMirror[1] - normalVec[1] * vectorBeforeMirror[0] #+ pointMirror1ZylKart[2]
  
  var alphaMirror = arcsin(abs(normalVec[0] * vectorBeforeMirror[0] + normalVec[1] * vectorBeforeMirror[1] + normalVec[2] * vectorBeforeMirror[2]) / (sqrt((normalVec[0] * normalVec[0] + normalVec[1] * normalVec[1] + normalVec[2] * normalVec[2])) * sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2]))))

  vectorAxis = vectorAxis / (sqrt((vectorAxis[0] * vectorAxis[0] + vectorAxis[1] * vectorAxis[1] + vectorAxis[2] * vectorAxis[2])))
  vectorBeforeMirror = vectorBeforeMirror / (sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2])))
  echo vectorAxis
  echo normalVec
  echo vectorBeforeMirror
  var vectorAfterMirror = vec3(0.0)
  vectorAfterMirror[0] = vectorBeforeMirror[0] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[1] * vectorAxis[2] - vectorBeforeMirror[2] * vectorAxis[1]) * sin( 2.0 * alphaMirror)
  vectorAfterMirror[1] = vectorBeforeMirror[1] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[2] * vectorAxis[0] - vectorBeforeMirror[0] * vectorAxis[2]) * sin( 2.0 * alphaMirror)
  vectorAfterMirror[2] = vectorBeforeMirror[2] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[0] * vectorAxis[1] - vectorBeforeMirror[1] * vectorAxis[0]) * sin( 2.0 * alphaMirror)

  var alphaTest = arccos((abs(vectorAfterMirror[0] * vectorBeforeMirror[0] + vectorAfterMirror[1] * vectorBeforeMirror[1] + vectorAfterMirror[2] * vectorBeforeMirror[2])) / (sqrt((vectorAfterMirror[0] * vectorAfterMirror[0] + vectorAfterMirror[1] * vectorAfterMirror[1] + vectorAfterMirror[2] * vectorAfterMirror[2])) * sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2]))))
  var pointAfterMirror = pointMirror + 200.0 * vectorAfterMirror
  var alphaVec = vec3(0.0)
  alphaVec[0] = radToDeg(alphaMirror)
  case pointOrAngle
  of "angle":
    result = alphaVec
  of "pointMirror":
    result = pointMirror
  of "pointAfter":
    result = pointAfterMirror
  of "vectorAfter":
    result = vectorAfterMirror]#

 


#[const
  RAYTRACER_LENGTH_COLDBORE = 9756.0
  RAYTRACER_LENGTH_PIPE_CB_VT3 = 2571.5
  RAYTRACER_LENGTH_PIPE_VT3_XRT = 150.0]#
var centerExitPipeVT3XRT = RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT
suite "Tests":
  test "Mirrors":
    let 
      d = 83.0 #mm ## distance between center of colbore at XRT and center of XRT (where the focal point is on the minus x axis)
      r1 = 83.405
      xSep = 4.284
      beta = degToRad(0.767)
      beta3 = 3.0 * beta
      
      lMirror = 225.0
      distanceMirrors = cos(beta3) * (xSep + lMirror)
      r2 = r1 - lMirror * sin(beta) #225mm is the length of the mirrors
      r3 = r2 - 0.5 * xSep * tan(beta)
      r4 = r3 - 0.5 * xSep * tan(beta3)
      r5 = r4 - lMirror * sin(beta3)
    var pointExitCB = vec3(0.0)
    pointExitCB[0] = - 1.5
    var pointEntranceXRT = vec3(0.0)
    pointEntranceXRT[2] = centerExitPipeVT3XRT
    pointEntranceXRT[0] = -1.5
    var pointEntranceXRTZylKart = vec3(0.0)
    pointEntranceXRTZylKart[0] = pointEntranceXRT[0] + d 
    pointEntranceXRTZylKart[1] = pointEntranceXRT[1]
    pointEntranceXRTZylKart[2] = pointEntranceXRT[2] - centerExitPipeVT3XRT 
  
    var pointExitCBZylKart = vec3(0.0)
    pointExitCBZylKart[0] = pointExitCB[0] + d 
    pointExitCBZylKart[1] = pointExitCB[1]
    pointExitCBZylKart[2] = pointExitCB[2] - centerExitPipeVT3XRT


    var pointMirror1 = findPosXRT(pointEntranceXRTZylKart, pointExitCBZylKart, r1, r2, beta, lMirror, 0.0,  0.001, 1.0, 1.1)
    
    
    var s = getVectoraAfterMirror(pointEntranceXRTZylKart, pointExitCBZylKart, pointMirror1, beta, "angle") ##poinExitPipeVT3 = pointEntranceXRT
    
    
    var vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRTZylKart, pointExitCBZylKart, pointMirror1, beta, "vectorAfter")
    var pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    var pointMirror2 = findPosXRT(pointAfterMirror1, pointMirror1, r4, r5, beta3, lMirror, distanceMirrors, 0.01, 0.0, 2.5)
    var t = getVectoraAfterMirror(pointAfterMirror1, pointMirror1, pointMirror2, beta3, "angle")
    check 2.0 * radToDeg(beta) - 0.5 * s[0] <= 0.5 * t[0] + 0.001 and 2.0 * radToDeg(beta) - 0.5 * s[0] >= 0.5 * t[0] - 0.001