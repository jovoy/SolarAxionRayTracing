import glm/vec
import math
# kompilieren und ausführen: nim cpp -r aEL.nim, nim c -r --threads:on --showAllMismatches:on aEL.nim # nim cpp -r --gc:boehm --verbosity:3 Raytracer2018.nim ##hdfview likelihood_2018_2.h5
# NOTE: Nein, einfach mit:
# nim c -d:danger -d:H5_FUTURE Raytracer2018.nim
import strutils
import algorithm
import plotly
import random
import sequtils, os
import nimhdf5
import chroma
#import ingrid/[tos_helpers, likelihood, ingrid_types]

import numericalnim, ggplotnim, strformat
import strscans
import streams, tables






##################rayTracer###############################

#degToRad(angle has to be done in Raytracer2014 for cos and sin



################################
# VARIABLES from rayTracer.h
const
  RAYTRACER_DISTANCE_SUN_EARTH =  1.5e14  #mm #ok
  radiusSun = 6.9e11 #mm #ok
  radiusCB = 21.5 #mm #ok
  RAYTRACER_LENGTH_COLDBORE = 9756.0 #mm half B field to end of CB #ok
  RAYTRACER_LENGTH_COLDBORE_9T = 9260.0 #mm half B field to half B field #ok
  RAYTRACER_LENGTH_PIPE_CB_VT3 = 2571.5 #mm should stay the same #from beam pipe drawings #ok
  radiusPipeCBVT3 = 39.64 #30.0 #mm smallest aperture between end of CB and VT4 # 43 mm but only 85% of area: 39,64mm #ok
  RAYTRACER_LENGTH_PIPE_VT3_XRT = 150.0 #mm from drawings #198.2 #mm from XRT drawing #ok
  radiusPipeVT3XRT = 35.0#25.0 #mm from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm #ok
  RAYTRACER_FOCAL_LENGTH_XRT = 1500.0#1485.0 #mm is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 #mm was the Telescope of 2014 (MPE XRT) also: Aperatur changed #ok
  distanceCBAxisXRTAxis = 0.0#62.1#58.44 #mm from XRT drawing #there is no difference in the axis even though the picture gets transfered 62,1mm down, but in the detector center
  RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW = 0.0 #mm #no change, because don't know # is actually -10.0 mm
  numberOfPointsEndOfCB = 200
  numberOfPointsSun = 20000 #100000
## Chipregions#####

const
  CHIPREGIONS_CHIP_X_MIN = 0.0
  CHIPREGIONS_CHIP_X_MAX = 14.0
  CHIPREGIONS_CHIP_Y_MIN = 0.0
  CHIPREGIONS_CHIP_Y_MAX = 14.0
  CHIPREGIONS_CHIP_CENTER_X = 7.0
  CHIPREGIONS_CHIP_CENTER_Y = 7.0
  CHIPREGIONS_GOLD_X_MIN = 4.5
  CHIPREGIONS_GOLD_X_MAX = 9.5
  CHIPREGIONS_GOLD_Y_MIN = 4.5
  CHIPREGIONS_GOLD_Y_MAX = 9.5
  CHIPREGIONS_SILVER_RADIUS_MAX = 4.5
  CHIPREGIONS_BRONZE_RADIUS_MAX = 5.5

################################



var fluxFractionGold = 0.0 #dies muss eine globale var sein
proc getFluxFractionGold(): float64 =

  result = fluxFractionGold

var fluxFractionSilver = 0.0
proc getFluxFractionSilver(): float64 =
  result = fluxFractionSilver

var fluxFractionBronze = 0.0
proc getFluxFractionBronze(): float64 =
  result = fluxFractionBronze

var fluxFractionDetector = 0.0
proc getFluxFractionDetector(): float64 =
  result = fluxFractionDetector

var fluxFractionTotal = 0.0
proc getFluxFractionTotal(): float64 =
  result = fluxFractionTotal


proc getFluxFraction(chipRegionstring:string): float64 =
  case chipRegionstring
  of "region: gold":
    result = getFluxFractionGold()
  of "region: silver":
    result = getFluxFractionSilver()
  of "region: bronze":
    result = getFluxFractionBronze()
  of "region: goldPlusSilver":
    result = getFluxFractionGold() + getFluxFractionSilver()
  of "region: goldPlusSilverPlusBronze":
    result = getFluxFractionGold() + getFluxFractionSilver() + getFluxFractionBronze()

  else: echo "Error: Unknown chip region!"

########################################### some functions we're gonna need later#################################################
# Now let's get some functions for a random point on a disk, a biased random point from the solar model and a biased random energy ##

proc getRandomPointOnDisk(center: Vec3, radius:float64) : Vec3 =

  ## This function gets a random point on a disk --> in this case this would be the entrance of the coldbore ##

  var
    x = 0.0
    y = 0.0
    r = radius * sqrt(random(1.0))  #random(1.0)
  # _randomGEnerator -> Circle(x,y,r)  is done through the following ### gives difference, since the root function Circle probably uses a differen random algorithm
    angle = 360 * random(1.0)  #random angle
  x = cos(degToRad(angle)) * r
  y = sin(degToRad(angle)) * r
  var vector = vec3(x,y,0.0)
  vector = vector + center
  result = vector

proc getEmRatesSums(energyVec: seq[float], emissionRates: seq[float]) : seq[float] =

  ## This function gives a vector, where thr first entry is the summ over all emissionrates (over all radii and energies) and the next are the for each radius over all energies plus the value of the last entry (except for the second entry) ##

  var
    sumAllEmRates = 0.0
    sumEmRatesR = 0.0
    emRateVecSums : seq[float]

  for i in 0..<emissionRates.len:
    sumAllEmRates = sumAllEmRates + emissionRates[i]
  emRateVecSums.add(sumAllEmRates)

  for r in 0..396:
    for iEnergy in 0..<energyVec.len:
      sumEmRatesR = sumEmRatesR + emissionRates[iEnergy + (r * energyVec.len)]
    emRateVecSums.add(sumEmRatesR)

  result = emRateVecSums


proc getRandomPointFromSolarModel(center : Vec3, radius : float64 , emRateVecSums : seq[float]) : Vec3 =

  ## This function gives the coordinates of a random point in the sun, biased by the emissionrates (which depend on the radius and the energy) ##

  var
    x = 0.0
    y = 0.0
    z = 0.0
    r :float
    angle1 = 360 * random(1.0)
    angle2 = 180 * random(1.0)
    i = random(emRateVecSums[0])

  for iRad in 1..<emRateVecSums.len:
    if iRad != 1 and i > emRateVecSums[iRad-1] and i <= emRateVecSums[iRad]:
      r =  (0.0015 + (iRad-1).float * 0.0005) * radius
    elif iRad == 1 and i >= 0.0 and i <= emRateVecSums[iRad]:
      r =  (0.0015 + (iRad-1).float * 0.0005) * radius


  x = cos(degToRad(angle1)) * sin(degToRad(angle2)) * r
  y = sin(degToRad(angle1)) * sin(degToRad(angle2)) * r
  z = cos(degToRad(angle2)) * r
  var vector = vec3(x,y,z)
  vector = vector + center

  result = vector

proc evenlyDistributedEnergies(energies :  seq[float], emissionRates : seq[float] , enOrEm: string) : seq[float] =

  ## get the energies and their emission rates as evenly distributed energies, because from the solar model we get very unevenly distributed energies and their emission rates ##

  var
    energyGl : float
    enrgiesGl: seq[float]
    emrateGl : float
    emratesGl : seq[float]
    x : float
    d_plus : float
    d_minus : float
  for iRad in 0..396:
    energyGl = 0.0
    emrateGl = 0.0
    for iEnergy in 0..<energies.len:
      energyGl = iEnergy.float * 0.0515 + 0.005
      enrgiesGl.add(energyGl)
      for iEnergy2 in 0..<energies.len:
        x = sqrt((energyGl-energies[iEnergy2])*(energyGl-energies[iEnergy2]))

        if iEnergy2 != 0 and iEnergy2 != (energies.len - 1):
          d_plus = sqrt((energyGl-energies[iEnergy2+1])*(energyGl-energies[iEnergy2+1]))
          d_minus = sqrt((energyGl-energies[iEnergy2-1])*(energyGl-energies[iEnergy2-1]))
          if x < d_minus and x < d_plus:
            if (energyGl-energies[iEnergy2]) < 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)] * (1- (x / sqrt((energies[iEnergy2]-energies[iEnergy2-1])*(energies[iEnergy2]-energies[iEnergy2-1])))) + emissionRates[(iEnergy2-1) + (iRad * energies.len)] * (1- (d_minus / sqrt((energies[iEnergy2]-energies[iEnergy2-1])*(energies[iEnergy2]-energies[iEnergy2-1]))))
              emratesGl.add(emrateGl)

            elif (energyGl-energies[iEnergy2]) > 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)] * (1- (x / sqrt((energies[iEnergy2]-energies[iEnergy2+1])*(energies[iEnergy2]-energies[iEnergy2+1])))) + emissionRates[(iEnergy2+1) + (iRad * energies.len)] * (1- (d_plus / sqrt((energies[iEnergy2]-energies[iEnergy2+1])*(energies[iEnergy2]-energies[iEnergy2+1]))))
              emratesGl.add(emrateGl)

        elif iEnergy2 == 0:
          d_plus = sqrt((energyGl-energies[iEnergy2+1])*(energyGl-energies[iEnergy2+1]))
          if x < d_plus:
            if (energyGl-energies[iEnergy2]) < 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)]  * (1- (x / sqrt((energies[iEnergy2]-0.0)*(energies[iEnergy2]-0.0)))) #+ (emissionRates[(iEnergy2+1) + (iRad * energies.len)] * sqrt((energyGl-energies[iEnergy2+1])*(energyGl-energies[iEnergy2+1])))
              emratesGl.add(emrateGl)

            elif(energyGl-energies[iEnergy2]) > 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)] * (1- (x / sqrt((energies[iEnergy2]-energies[iEnergy2+1])*(energies[iEnergy2]-energies[iEnergy2+1])))) + emissionRates[(iEnergy2+1) + (iRad * energies.len)] * (1- (d_plus / sqrt((energies[iEnergy2]-energies[iEnergy2+1])*(energies[iEnergy2]-energies[iEnergy2+1]))))
              emratesGl.add(emrateGl)
        elif iEnergy2 == (energies.len - 1):
          d_minus = sqrt((energyGl-energies[iEnergy2-1])*(energyGl-energies[iEnergy2-1]))
          if x < d_minus:
            if (energyGl-energies[iEnergy2]) < 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)] * (1- (x / sqrt((energies[iEnergy2]-energies[iEnergy2-1])*(energies[iEnergy2]-energies[iEnergy2-1])))) + emissionRates[(iEnergy2-1) + (iRad * energies.len)] * (1- (d_minus / sqrt((energies[iEnergy2]-energies[iEnergy2-1])*(energies[iEnergy2]-energies[iEnergy2-1]))))
              emratesGl.add(emrateGl)
            elif (energyGl-energies[iEnergy2]) > 0.0:
              emrateGl = emissionRates[(iEnergy2) + (iRad * energies.len)] * (1- (x / sqrt((energies[iEnergy2]-0.0)*(energies[iEnergy2]-0.0)))) #+ (emissionRates[(iEnergy2-1) + (iRad * energies.len)] * (energyGl-energies[iEnergy2-1]))
              emratesGl.add(emrateGl)
  if enOrEm == "energies":
    result = enrgiesGl
  elif enOrEm == "emrates":
    result = emratesGl

proc getRandomEnergyFromSolarModel(emRateVecSums : seq[float], vectorInSun: Vec3, center : Vec3, radius : float64, energies :  seq[float],emissionRates: seq[float]) : float =

  ## This function gives a random energy for an event at a given radius, biased by the emissionrates at that radius. This only works if the energies to choose from are evenly distributed ##

  var
    rad = sqrt(vectorInSun[0]*vectorInSun[0]+vectorInSun[1]*vectorInSun[1]+ (vectorInSun[2]-center[2])*(vectorInSun[2]-center[2]))
    r = rad / radius
    iRad : int
    energy:float
    emRateEnergySum : float
    emRateEnergySumPrev : float
    i : float
    emRateEnergySumAll : float

  var
    indexRad = (r - 0.0015) / 0.0005

  if indexRad - 0.5 > floor(indexRad):
    iRad = int(ceil(indexRad))
  else: iRad = int(floor(indexRad))

  for iE in 0..<233:
    emRateEnergySumAll = emRateEnergySumAll + emissionRates[(iE) + (iRad * 233)]
  #echo emRateEnergySumAll


  if iRad != 0:
    #i = random((emRateVecSums[iRad+1] - emRateVecSums[iRad])) # no because this is not with energy Gl
    i = random(emRateEnergySumAll)
  else:
    i = random(emRateVecSums[iRad+1])

  for iEnergy in 0..<233: #energise.len/ radii
    emRateEnergySumPrev = emRateEnergySum
    emRateEnergySum = emRateEnergySumPrev + emissionRates[(iEnergy) + (iRad * 233)]#energise.len/ radii
    if i > emRateEnergySumPrev and i <= emRateEnergySum:
      energy = energies[iEnergy]
    #if iEnergy == 232:
      #echo "sum over all energies 2"
      #echo emRateEnergySum

  result = energy

## The following are some functions to determine inetersection of the rays with the geometry (pipes, magnet, mirrors) of the setup ##

proc lineIntersectsCircle(point_1 : Vec3, point_2 : Vec3, center : Vec3, radius : float64, intersect : Vec3) : bool = # probably still some error (lambda1 -> infinity)

  ## Now a function to see if the linesfrom the sun will actually intersect the circle area from the magnet entrance (called coldbore) etc. ##

  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  var intersect = vec3(0.0)
  intersect = point_1 + lambda1 * vector
  var r_xy_intersect = sqrt(intersect[0] * intersect[0] + intersect[1] * intersect[1])
  if r_xy_intersect < radius:
    return true
  else:
    return false



proc lineIntersectsCylinderOnce(point_1 : Vec3, point_2 : Vec3, centerBegin : Vec3, centerEnd : Vec3, radius : float64, intersect : Vec3) : bool =

  ## Also a function to know if the line intersected at least the whole magnet, and then only once, because else the axions would have just flown through ##

  var vector = vec3(0.0)
  var intersect = vec3(0.0)
  vector = point_2 - point_1
  var lambda_dummy : float64
  lambda_dummy = ( -1000.0 - point_1[2] ) / vector[2]
  var dummy = vec3(0.0)
  dummy = point_1 + lambda_dummy * vector
  var vector_dummy = vec3(0.0)
  vector_dummy = point_2 - dummy
  var factor : float64
  factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
  var p : float64
  p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
  var q : float64
  q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - radius*radius) / factor
  var lambda_1 : float64
  var lambda_2 : float64
  lambda_1 = -p/2.0 + sqrt( p*p/4.0 - q)
  lambda_2 = -p/2.0 - sqrt( p*p/4.0 - q)
  var intersect_1 = vec3(0.0)
  intersect_1 = dummy + lambda_1 * vector_dummy
  var intersect_2 = vec3(0.0)
  intersect_2 = dummy + lambda_2 * vector_dummy
  var intersect_1_valid : bool
  var intersect_2_valid : bool
  intersect_1_valid = (intersect_1[2] > centerBegin[2] ) and (intersect_1[2] < centerEnd[2])
  intersect_2_valid = (intersect_2[2] > centerBegin[2] ) and (intersect_2[2] < centerEnd[2])
  if ( (intersect_1_valid and intersect_2_valid) or (not intersect_1_valid and not intersect_2_valid) ):
    return false
  elif (intersect_1_valid):
    intersect = intersect_1
    return true
  else:
    intersect = intersect_2
    return true

proc getIntersectLineIntersectsCylinderOnce(point_1 : Vec3, point_2 : Vec3, centerBegin : Vec3, centerEnd : Vec3, radius : float64, intersect : Vec3) : Vec3 =
  var vector = vec3(0.0)
  var intersect = vec3(0.0)
  vector = point_2 - point_1
  var lambda_dummy : float64
  lambda_dummy = ( -1000.0 - point_1[2] ) / vector[2]
  var dummy = vec3(0.0)
  dummy = point_1 + lambda_dummy * vector
  var vector_dummy = vec3(0.0)
  vector_dummy = point_2 - dummy
  var factor : float64
  factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
  var p : float64
  p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
  var q : float64
  q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - radius*radius) / factor
  var lambda_1 : float64
  var lambda_2 : float64
  lambda_1 = -p/2.0 + sqrt( p*p/4.0 - q)
  lambda_2 = -p/2.0 - sqrt( p*p/4.0 - q)
  var intersect_1 = vec3(0.0)
  intersect_1 = dummy + lambda_1 * vector_dummy
  var intersect_2 = vec3(0.0)
  intersect_2 = dummy + lambda_2 * vector_dummy
  var intersect_1_valid : bool
  var intersect_2_valid : bool
  intersect_1_valid = (intersect_1[2] > centerBegin[2] ) and (intersect_1[2] < centerEnd[2])
  intersect_2_valid = (intersect_2[2] > centerBegin[2] ) and (intersect_2[2] < centerEnd[2])
  if (intersect_1_valid):
    intersect = intersect_1
  else:
    intersect = intersect_2
  result = intersect



proc circleEdges( theR1 : seq[float]): seq[seq[float]] =

  ## A function to create the coordinates of the mirrors of the telescope (they're made of glass) to substract from the overall picture, because Xrays in that range don't pass through glass ##

  const
    sizeViewfield = 48.0 #mm
    d = 83.0 #mm

  var
    circleEdge = vec3(0.0)
    circleEdges = newSeqWith(1401, newSeq[float64](1401))

  for R1 in theR1:
    for phi in 34500 ..< 37500:
      for dgl in 0 .. 20:
        circleEdge[1]= (R1 + (float(dgl)*0.01)) * sin(degToRad(float(phi)*0.01))
        circleEdge[0]= (R1 + (float(dgl)*0.01)) * cos(degToRad(float(phi)*0.01)) - d
        var coord_X = floor(circleEdge[0] / (sizeViewfield/1400.0)) + 700
        var coord_Y = floor(circleEdge[1] / (sizeViewfield/1400.0)) + 700
        if coord_X >= 0.0:
          if coord_Y >= 0.0:
            if coord_X <= 1400.0:
              if coord_Y <= 1400.0:
                circleEdges[int(coord_X)][int(coord_Y)] = circleEdges[int(coord_X)][int(coord_Y)] + 1

  result = circleEdges

proc getPixelValue(intersects : Vec3): Vec3 =
  const sizeViewfield = 48.0 #mm
  var intersectsPix = vec3(0.0)
  intersectsPix[0] = floor(intersects[0] / (sizeViewfield/1400.0)) + 700
  intersectsPix[1] = floor(intersects[1] / (sizeViewfield/1400.0)) + 700
  result = intersectsPix


proc lineIntersectsCircleEdge(circleEdges : seq[seq[float64]]  ,  intersectsPix: Vec3): bool =
  var
    coord_X = int(intersectsPix[0])
    coord_Y = int(intersectsPix[1])
  if circleEdges[coord_X][coord_Y] == 0.0 :
    return false
  else: return true

proc lineIntersectsArea( R1 : float64, prevR1 : float64, intersect : Vec3): bool =
  const
    d = 83.0 #mm
  var r_xy_intersect = sqrt(intersect[1] * intersect[1] + (intersect[0]+d) * (intersect[0]+d))
  if r_xy_intersect <= R1 and r_xy_intersect >= (prevR1) :
    return true
  else:
    return false

## Some functions to include files from outside like the run file and the emissionrate/energy files ##

proc getXandY( h5file : string , dsetgrp1 :string, numFirstRun : int, numLastRun : int, chip : string, xOrY : string) : seq[float] =
  ## get the x and y values from the run-file to compare them to our model ##
  echo fileExists(h5file)
  var h5f = H5file(h5file, "r")
  for grp in items(h5f, start_path = dsetgrp1):
    if grp.name / chip in h5f:
      echo grp
      #var energy = h5f[(grp.name / chip / "energyFromCharge"), float64]
      if xOrY == "X":
        result.add h5f[(grp.name / chip / "centerX"), float64]
      else:
        result.add h5f[(grp.name / chip / "centerY"), float64]
  #if xOrY == "X":
  #  result = valuesX
  #elif xOrY == "Y":
  #  result = valuesY
  #else: return @[0.0]

proc getAxionEnergy(): seq[float] =
  var
    f = open("energies.txt")
    line = ""
    energy : seq[float]
    energystring : string

  #doAssert f.getPosition()
  while f.readLine(line):
    let energystring = line
    energy.add(parseFloat(energystring))

  defer: f.close()
  return energy

proc getfluxfracVec () : seq[float] =
  var
    f = open("flux_fractions_with_fB.txt") # with Bose-Einstein distribution of the thermal photon bath
    line = ""
    fluxstring : string
    fluxfractions : seq[float]

  while f.readLine(line):
    let fluxstring = line
    fluxfractions.add(parseFloat(fluxstring))

  defer: f.close()
  return fluxfractions

proc getfluxfrac (radius : float, energyVec: seq[float], iEnergy : int) : int =
  var
    linenum : int
    linenumfloat : float
    posfluxfraction : int
    radiusPerc : float

  radiusPerc = radius / radiusSun
  linenumfloat = ((radiusPerc - 0.0015) / 0.0005) + 1.0

  if round(linenumfloat, 1) < (floor(linenumfloat) + 0.5) :
    linenum = int(floor(linenumfloat))
  else:
    linenum = int(ceil(linenumfloat))

  if linenum > 0 and linenum < 398:
    posfluxfraction = (linenum - 1 ) * (energyVec.len ) + (iEnergy )
  elif linenum > 397:
    posfluxfraction = (linenum - 2 ) * (energyVec.len ) + (iEnergy )
  elif linenum == 0:
    posfluxfraction = (linenum ) * (energyVec.len ) + (iEnergy )
  elif linenum < 0:
    posfluxfraction = (linenum + 2) * (energyVec.len ) + (iEnergy )

  return posfluxfraction#linenum##posfluxfraction#

proc getEmRateVec() : seq[float] =
  var
    f = open("emission_rates_Hz.txt")
    line = ""
    emissionstring : string
    emissionrates : seq[float]

  while f.readLine(line):
    let emissionstring = line
    emissionrates.add(parseFloat(emissionstring))

  defer: f.close()
  return emissionrates


proc findPosXRT*(pointXRT : Vec3, pointCB : Vec3, r1 : float, r2 : float, angle : float, lMirror : float, distMirr : float, uncer : float, sMin : float, sMax : float) : Vec3 =

  ##this is to find the position the ray hits the mirror shell of r1. it is after transforming the ray into a cylindrical coordinate system, that has the middle and the beginning of the mirror "cylinders" as its origin and is in cylindrical coordinates
    
  var
    point = pointCB
    s : float
    term : float
    sMinHigh = (sMin * 100000.0).int
    sMaxHigh = (sMax * 100000.0).int
    pointMirror = vec3(0.0)
  let direc = pointXRT - pointCB
    
      
  for i in sMinHigh..sMaxHigh:
    s = i.float / 100000.0
    term = sqrt((point[0] + s * direc[0]) * (point[0] + s * direc[0]) + (point[1] + s * direc[1]) *  (point[1] + s * direc[1])) - ((r2 - r1) * (point[2] + s * direc[2] - distMirr) / (cos(angle) * lMirror))
        
    if r1 < term + uncer and r1 > term - uncer:
      pointMirror = point + s * direc ## sometimes there are 3 different s for which this works, in that case the one with the highest value is taken
      result = pointMirror


proc getVectoraAfterMirror*(pointXRT : Vec3, pointCB : Vec3, pointMirror : Vec3, angle : float, pointOrAngle : string) : Vec3 =

  ## this is to find the vector after the reflection on the respective mirror
      
  var normalVec = vec3(0.0)
  normalVec[0] = pointMirror[0]
  normalVec[1] = pointMirror[1]
  normalVec[2] = tan(angle) * sqrt(pointMirror[0] * pointMirror[0] + pointMirror[1] * pointMirror[1])
  var vectorBeforeMirror =  pointXRT - pointCB
  var vectorAxis = vec3(0.0) ## this is the vector product of the normal vector on pointMirror pointing in the direction of the radius of the cylinder and the vector of the ray
  vectorAxis[0] = normalVec[1] * vectorBeforeMirror[2] - normalVec[2] * vectorBeforeMirror[1] 
  vectorAxis[1] = normalVec[2] * vectorBeforeMirror[0] - normalVec[0] * vectorBeforeMirror[2] 
  vectorAxis[2] = normalVec[0] * vectorBeforeMirror[1] - normalVec[1] * vectorBeforeMirror[0] 
      
  var alphaMirror = arcsin(abs(normalVec[0] * vectorBeforeMirror[0] + normalVec[1] * vectorBeforeMirror[1] + normalVec[2] * vectorBeforeMirror[2]) / (sqrt((normalVec[0] * normalVec[0] + normalVec[1] * normalVec[1] + normalVec[2] * normalVec[2])) * sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2]))))
      
  vectorAxis = vectorAxis / (sqrt((vectorAxis[0] * vectorAxis[0] + vectorAxis[1] * vectorAxis[1] + vectorAxis[2] * vectorAxis[2])))
  vectorBeforeMirror = vectorBeforeMirror / (sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2])))
  var vectorAfterMirror = vec3(0.0)
  vectorAfterMirror[0] = vectorBeforeMirror[0] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[1] * vectorAxis[2] - vectorBeforeMirror[2] * vectorAxis[1]) * sin( 2.0 * alphaMirror)
  vectorAfterMirror[1] = vectorBeforeMirror[1] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[2] * vectorAxis[0] - vectorBeforeMirror[0] * vectorAxis[2]) * sin( 2.0 * alphaMirror)
  vectorAfterMirror[2] = vectorBeforeMirror[2] * cos(2.0 * alphaMirror) - (vectorBeforeMirror[0] * vectorAxis[1] - vectorBeforeMirror[1] * vectorAxis[0]) * sin( 2.0 * alphaMirror)
      
  var alphaTest = arccos((abs(vectorAfterMirror[0] * vectorBeforeMirror[0] + vectorAfterMirror[1] * vectorBeforeMirror[1] + vectorAfterMirror[2] * vectorBeforeMirror[2])) / (sqrt((vectorAfterMirror[0] * vectorAfterMirror[0] + vectorAfterMirror[1] * vectorAfterMirror[1] + vectorAfterMirror[2] * vectorAfterMirror[2])) * sqrt((vectorBeforeMirror[0] * vectorBeforeMirror[0] + vectorBeforeMirror[1] * vectorBeforeMirror[1] + vectorBeforeMirror[2] * vectorBeforeMirror[2]))))
  var pointAfterMirror = pointMirror + 200.0 * vectorAfterMirror
  var alphaVec = vec3(0.0)
  alphaVec[0] = radToDeg(alphaTest)
  case pointOrAngle
  of "angle":
    result = alphaVec
  of "pointMirror":
    result = pointMirror
  of "pointAfter":
    result = pointAfterMirror
  of "vectorAfter":
    result = vectorAfterMirror


## Now some functions for the graphs later, that store the data in heatmaps and then give them out with plotly ##


proc prepareheatmap(numberofrows : int, numberofcolumns : int,
                    start_x : float, stop_x : float, start_y : float,
                    stop_y : float,
                    data_X : seq[float], data_Y : seq[float], weight1 : seq[float],
                    norm : float64): seq[seq[float]] =
  ## This function prepares a heatmap out of given X and Y values with the z value (the number of entries in a certain pixel) as the weight of the event of the X and Y value ##

  var stepsize_X = 0.0  # number of colums is the number of entries in the array in the seq and number of rows is the number of arrays in the seq
  stepsize_X = (stop_x - start_x)/float(numberofrows)
  var stepsize_Y = 0.0
  stepsize_Y = (stop_y - start_y)/float(numberofcolumns)
  #[
    |-------|-------|-------|
    |-------|-------|-------|
    |-------|-------|-------|
    |-------|-------|-------|
    |-------|-------|-------|
  ]#
  var heatmaptable = newSeqWith(numberofrows, newSeq[float](numberofcolumns))
  for i, value in data_X:
    var coord_X = floor((data_X[i] - start_x) / stepsize_X)
    var coord_Y = floor((data_Y[i] - start_y) / stepsize_Y)
    if coord_X >= 0:
      if coord_Y >= 0:
        if coord_X <= float(numberofrows):
          if coord_Y <= float(numberofcolumns):
            heatmaptable[int(coord_Y)][int(coord_X)] = heatmaptable[int(coord_Y)][int(coord_X)] + 1*weight1[i]/norm
  result = heatmaptable

proc getMaxVal(table : seq[seq[float]], numberofrows : int) : float =
  var maxVals : seq[float]
  var maxVal : float64
  for i in 0 ..< numberofrows:
    maxVals.add(max(table[i]))
  maxVal = max(maxVals)
  result = maxVal

proc getLenght(table : seq[seq[float]], numberofrows : int) : int =
  var lengths : seq[float]
  var length : int
  for i in 0 ..< numberofrows:
    lengths.add(max(table[i]))
  length = lengths.len
  result = length

proc drawfancydiagrams(diagramtitle : string, objectstodraw : seq[seq[float]], width : int) : float =

  ## this function draws a hdiagram out a given heatmap ##

  let
    # The GL heatmap is also supported as HeatMapGL
    d = Trace[float32](mode: PlotMode.Lines, `type`: PlotType.HeatMap)

  d.colormap = ColorMap.Viridis
  d.zs = newSeqWith(width, newSeq[float32](width))
  for x in 0 ..< width:
    for y in 0 ..< width:
      if x < width:
        if y < width:
          if x > 0:
            if y > 0:
              d.zs[y][x] = objectstodraw[y][x]

  const
    y = @[float32(CHIPREGIONS_GOLD_Y_MIN * 3000.0 / 14.0), float32(CHIPREGIONS_GOLD_Y_MIN * 3000.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 3000.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 3000.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MIN * 3000.0 / 14.0)]
    x = @[float32(CHIPREGIONS_GOLD_X_MIN * 3000.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 3000.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 3000.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 3000.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 3000.0 / 14.0)]
  let
    d4 = Trace[float32](mode: PlotMode.LinesMarkers, `type`: PlotType.ScatterGL, ys: y, xs : x)


  let
    layout = Layout(title: diagramtitle, width: 800, height: 800,
                    xaxis: Axis(title: "x-axis [mm]"),#,  range: (0.0, 14.0)),
                    yaxis: Axis(title: "y-axis [mm]"), autosize: false)
    p = Plot[float32](layout: layout, traces: @[d, d4])
  echo p.save()
  p.show()


proc drawgraph(diagramtitle : string, data_X : seq[float], data_Y: seq[float], energyOrRadius: string) : float =

  ## This function draws a graph out of given x and y values ##

  var
    #colors : seq[Color]
    #colors = new_seq[Color](data_X.len)
    size = @[16.0]
    color = @[Color(r: 0.9, g: 0.1, b: 0.1, a: 1.0)]
    layout = Layout()
    #sizes = new_seq[float64](data_X.len)
  #for i in 0..<data_X.len:
    #sizes[i] = (10.0)
    #colors[i] = Color(r: 0.9, g: 0.1, b: 0.1, a: 1.0)#0xFAF0E6)#

  let d = Trace[float64](mode: PlotMode.Markers, `type`: PlotType.ScatterGL,
                             xs: data_X, ys: data_Y)
  #d.marker = Marker[float64](size: size, color: color)
  if energyOrRadius == "energy" :
    layout = Layout(title: diagramtitle, width: 800, height: 800,
                      xaxis: Axis(title: "Energy [keV]"),
                      yaxis: Axis(title: "Emission Rate"), autosize: false) #Fluxfraction [10^20 m^-2 year^-1 keV^-1]

  elif energyOrRadius == "radius" :
    layout = Layout(title: diagramtitle, width: 800, height: 800,
                      xaxis: Axis(title: "Radius [% of the radius of the sun]"),
                      yaxis: Axis(title: "Fluxfraction [10^20 m^-2 year^-1 keV^-1]"), autosize: false)
  let p = Plot[float](layout: layout, traces: @[d])
  echo p.save()
  p.show()


############done with the functions, let's use them############

## First lets access the solar model and calculate some necessary values
const solarModel = "./ReadSolarModel/resources/AGSS09_solar_model_stripped.dat"

var df = readSolarModelDf(solarModel)
df = df.filter(f{"Radius" <= 0.2})
echo df.pretty(precision = 10)

# to read a single column, e.g. radius:
#echo df1["Rho"].len

# now let's plot radius against temperature colored by density
ggplot(df, aes("Radius", "Temp", color = "Rho")) +
  geom_line() +
  ggtitle("Radius versus temperature of solar mode, colored by density") +
  ggsave("radius_temp_density.pdf")

var n_Z = newSeqWith(df["Rho"].len, newSeq[float](29)) #29 elements
var n_e : seq[float]
var distTemp : float
var temperature : int
var temperatures : seq[int]
let atomicMass = [1.0078,4.0026,3.0160,12.0000,13.0033,14.0030,15.0001,15.9949,16.9991,17.9991,20.1797,22.9897,24.3055,26.9815,28.085,30.9737,32.0675,35.4515,39.8775,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.9380,55.845,58.9331,58.6934] #all the 29 elements from the solar model file
let elements = ["H1", "He4","He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
# var solarTable = readSolarModel(solarModel)
let amu = 1.6605e-24 #grams
echo df["Mass"][6]
for iRadius in 0..< df["Rho"].len:
  n_Z[iRadius][1] = (df[elements[0]][iRadius].toFloat / atomicMass[0]) * (df["Rho"][iRadius].toFloat / amu) # Hydrogen
  for iZmult in 1..3:
    n_Z[iRadius][iZmult * 2] = ((df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat) / ((atomicMass[iZmult * 2 - 1] * df[elements[iZmult * 2 - 1]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)) + (atomicMass[iZmult * 2] * df[elements[iZmult * 2]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
    n_Z[iRadius][8] = ((df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat) / ((atomicMass[7] * df[elements[7]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[8] * df[elements[8]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[9] * df[elements[9]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
  for iZ in 10..<29:
    n_Z[iRadius][iZ] = (df[elements[iZ]][iRadius].toFloat / atomicMass[iZ]) * (df["Rho"][iRadius].toFloat / amu) # The rest
  n_e.add((df["Rho"][iRadius].toFloat/amu) * (1 + df[elements[0]][iRadius].toFloat/2))
  #echo log(parseFloat(solarTable["Temp"][iRadius]), 10.0) / 0.025
  for iTemp in 0..90:
    distTemp = log(df["Temp"][iRadius].toFloat, 10.0) / 0.025 - float(140 + 2 * iTemp)
    if abs(distTemp) <= 1.0: 
      temperature = 140 + 2 * iTemp
      #echo distTemp
  temperatures.add(temperature)
echo temperatures



echo n_e
echo n_Z



proc calculateFluxFractions(axionRadiationCharacteristic: string,
  detectorWindowAperture :float64 = 14.0,
  sunMisalignmentH :float64 = 0.0, #given in arcmin(?)
  sunMisalignmentV:float64 = 0.0,
  detectorMisalignmentX:float64 = 0.0,
  detectorMisalignmentY:float64 = 0.0,
  coldboreBlockedLength:float64 = 0.0) : int =

  var misalignmentSun = vec3(0.0) #works, but maybe bad package?
  misalignmentSun[2] = 0
  misalignmentSun[0] = tan(degToRad(sunMisalignmentH / 60)) * RAYTRACER_DISTANCE_SUN_EARTH
  misalignmentSun[1] = tan(degToRad(sunMisalignmentV / 60)) * RAYTRACER_DISTANCE_SUN_EARTH

  var misalignmentDetector = vec3(0.0)
  misalignmentDetector[0] = detectorMisalignmentX
  misalignmentDetector[1] = detectorMisalignmentY
  misalignmentDetector[2] = 0

  var centerSun = vec3(0.0)
  centerSun[0] = 0
  centerSun[1] = 0
  centerSun[2] = RAYTRACER_DISTANCE_SUN_EARTH

  centerSun = centerSun + misalignmentSun

  var radiusSunTwentyPecent = radiusSun * 0.2

  var centerEntranceCB = vec3(0.0)
  centerEntranceCB[0] = 0
  centerEntranceCB[1] = 0
  centerEntranceCB[2] = coldboreBlockedLength

  var centerExitCBMagneticField = vec3(0.0)
  centerExitCBMagneticField[0] = 0
  centerExitCBMagneticField[1] = 0
  centerExitCBMagneticField[2] = RAYTRACER_LENGTH_COLDBORE_9T

  var centerExitCB = vec3(0.0)
  centerExitCB[0] = 0
  centerExitCB[1] = 0
  centerExitCB[2] = RAYTRACER_LENGTH_COLDBORE

  var centerExitPipeCBVT3 = vec3(0.0)
  centerExitPipeCBVT3[0] = 0
  centerExitPipeCBVT3[1] = 0
  centerExitPipeCBVT3[2] = RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3

  var centerExitPipeVT3XRT = vec3(0.0)
  centerExitPipeVT3XRT[0] = 0
  centerExitPipeVT3XRT[1] = 0
  centerExitPipeVT3XRT[2] = RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT

  var
    integralNormalisation = 0.0
    integralTotal = 0.0
    integralDetector = 0.0
    integralBronze = 0.0
    integralSilver = 0.0
    integralGold = 0.0

  var misalignment = (sunMisalignmentH != 0.0) or (sunMisalignmentV != 0.0) or (detectorMisalignmentX != 0.0) or (detectorMisalignmentY != 0.0)




  var
    pointdataX : seq[float]
    pointdataY : seq[float]
    pointdataXBefore : seq[float]
    pointdataYBefore : seq[float]
    weights : seq[float]
    pixvalsX : seq[float]
    pixvalsY : seq[float]
  const
    allR1 = @[60.7095, 63.006 , 65.606 , 68.305 , 71.105 , 74.011 , 77.027 , 80.157 , 83.405 , 86.775 , 90.272 , 93.902 , 97.668 , 101.576 , 105.632]  ## the radii of the shells
    allXsep = @[4.0, 4.171, 4.140, 4.221, 4.190, 4.228, 4.245, 4.288, 4.284, 4.306, 4.324, 4.373, 4.387, 4.403, 4.481]
    #allR2 = @[0.0, 60.731, 63.237, 65.838, 68.538, 71.339, 74.246, 77.263, 80.394, 83.642]
    allAngles = @[0.0, 0.579, 0.603, 0.628, 0.654, 0.680, 0.708, 0.737, 0.767, 0.798, 0.830, 0.863, 0.898, 0.933, 0.970] ## the angles of the mirror shells coresponding to the radii above
    lMirror = 225.0 #mm Mirror length
    circleX = circleEdges( allR1)
    circleY = circleEdges( allR1)
    circleTotal = circleEdges( allR1)
    d = 83.0 #mm ## distance between center of colbore at XRT and center of XRT (where the focal point is on the minus x axis)
  let
    energy = getAxionEnergy()
    fluxfracs = getfluxfracVec()
    emrates = getEmRateVec()
    emratesSums = getEmRatesSums(energy, emrates)

  var
    fluxfractionrad : seq[float]
    fluxfractionen : seq[float]
    radii : seq[float]
    energies : seq[float]
    fluxfracints : seq[float]

  let
    energiesGl = evenlyDistributedEnergies(energy, emrates ,"energies")
    emratesGl = evenlyDistributedEnergies(energy, emrates ,"emrates")
  var ffff : seq[float]

  #for iRad in 0..396:
    #var
      #r= 0.0015 + (iRad.float * 0.0005)

      #fluxfr = fluxfracs[getfluxfrac(r, energy, iEnergy)]

    #fluxfracint = fluxfracint + (fluxfr * r)
    #var vecSun = getRandomPointFromSolarModel(centerSun,radiusSun,emratesSums)
    #radii.add(getRandomEnergyFromSolarModel(emratesSums, vecSun,centerSun, radiusSun,energy, emrates))
    #fluxfractionrad.add(1.0 + (iRad.float * 0.001))
  # count through all points in the magnet end (coldbore) and then in the sun (set to 1000 points), to connect them and get 1000000 lines (axion traces)


    #fluxfracints.add(fluxfracint)
    #energies.add(energyAx)
  for iRad in 0..2:
    var 
      r = 396.0 * iRad.float / 2.0
      i = r.int
    for iEnergy in 0..<energy.len:
      energies.add(energy[iEnergy])
      fluxfractionen.add(emrates[(iEnergy) + (i * 233)])








  for iSun in 1..numberOfPointsSun:
    integralNormalisation = integralNormalisation + 1
    var pointInSun = vec3(0.0)
    case  axionRadiationCharacteristic
    of "axionRadiation::characteristic::sar":
      pointInSun = getRandomPointFromSolarModel(centerSun,radiusSun,emratesSums)#getRandomPointOnDisk(centerSun,radiusSunTwentyPecent)#pointInSun = getRandomPointFromSolarModel(centerSun,radiusSun)
    of "axionRadiation::characteristic::def":
      echo "Error: Default radiation characteristic not implemented"
    else:
      echo "Error: Unknown axion radiation characteristic"

    var pointExitCBMagneticField = vec3(0.0)
    pointExitCBMagneticField = getRandomPointOnDisk(centerExitCBMagneticField, radiusCB)
    var radiusInSun = sqrt(pointInSun[0] * pointInSun[0] + pointInSun[1] * pointInSun[1])


    var vecSun = getRandomPointFromSolarModel(centerSun,radiusSun,emratesSums)
    var energyAx = getRandomEnergyFromSolarModel(emratesSums, vecSun,centerSun, radiusSun,energiesGl, emratesGl)
    var intersect = vec3(0.0)
    var pathCB : float64
    var intersectsEntranceCB : bool
    intersectsEntranceCB = lineIntersectsCircle(pointInSun, pointExitCBMagneticField, centerEntranceCB, radiusCB, intersect)
    var intersectsCB = false
    #echo energyAx

    if (not intersectsEntranceCB):
      intersectsCB = lineIntersectsCylinderOnce(pointInSun,pointExitCBMagneticField,centerEntranceCB,centerExitCBMagneticField,radiusCB,intersect)

    if (not intersectsEntranceCB and not intersectsCB): continue

    if (not intersectsEntranceCB): #generates problems with the weight because the weight is multiplied with the difference of the leght of the path of the particle and the legth of the coldbore
      intersect = getIntersectLineIntersectsCylinderOnce(pointInSun,pointExitCBMagneticField,centerEntranceCB,centerExitCBMagneticField,radiusCB,intersect)#pointInSun + ((centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
    #if (not intersectsCB):
      #intersect = pointInSun + ((centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
    pathCB = pointExitCBMagneticField[2] - intersect[2]
    var pointExitCB = vec3(0.0)

    if (not lineIntersectsCircle(pointInSun,pointExitCBMagneticField,centerExitCB,radiusCB,pointExitCB)): continue

    pointExitCB = pointInSun + ((centerExitCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)

    var pointExitPipeCBVT3 =vec3(0.0)

    if (not lineIntersectsCircle(pointExitCBMagneticField,pointExitCB,centerExitPipeCBVT3,radiusPipeCBVT3,pointExitPipeCBVT3)) : continue

    pointExitPipeCBVT3 = pointExitCBMagneticField + ((centerExitPipeCBVT3[2] - pointExitCBMagneticField[2]) / (pointExitCB - pointExitCBMagneticField)[2]) * (pointExitCB - pointExitCBMagneticField)
    var pointExitPipeVT3XRT = vec3(0.0)#seq[float]

    if (not lineIntersectsCircle(pointExitCB,pointExitPipeCBVT3,centerExitPipeVT3XRT,radiusPipeVT3XRT,pointExitPipeVT3XRT)) : continue

    pointExitPipeVT3XRT = pointExitCB + ((centerExitPipeVT3XRT[2] - pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) * (pointExitPipeCBVT3 - pointExitCB)

    #pointExitPipeVT3XRT = pointInSun - pointExitCBMagneticField
    #pointExitPipeVT3XRT = pointExitPipeVT3XRT/ sqrt(pointExitPipeVT3XRT[0]*pointExitPipeVT3XRT[0] + pointExitPipeVT3XRT[1]*pointExitPipeVT3XRT[1] + pointExitPipeVT3XRT[2]*pointExitPipeVT3XRT[2])
    #pointExitPipeVT3XRT = pointExitCBMagneticField + pointExitPipeVT3XRT * (RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT)

    var vectorBeforeXRT = vec3(0.0)
    vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB




    ###################von CB (coldbore(pipe in Magnet)) zum XRT (XrayTelescope)#######################

    var pointEntranceXRT = vec3(0.0)
    pointEntranceXRT[0] = pointExitPipeVT3XRT[0] #- distanceCBAxisXRTAxis
    pointEntranceXRT[1] = pointExitPipeVT3XRT[1]
    pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

    if (getPixelValue(pointEntranceXRT)[0] > 1400.0 or getPixelValue(pointEntranceXRT)[1] > 1400.0): continue
    if lineIntersectsCircleEdge(circleTotal, getPixelValue(pointEntranceXRT)): continue

    #echo getPixelValue(pointEntranceXRT)
    #echo lineIntersectsCircleEdge(circleTotal, getPixelValue(pointEntranceXRT))
    # x and y is interchanged from this point on since it somehow is interchanged in the picture
    if  pointEntranceXRT[1] <= 1.0 and pointEntranceXRT[1] >= -1.0: continue #there is a 2mm wide graphit block between each glass mirror, to seperate them in the middle of the x-Ray telescope

    var
      vectorAfterXRTCircular = vec3(0.0)
      radius1 = sqrt((pointEntranceXRT[0]+d) * (pointEntranceXRT[0]+d) + (pointEntranceXRT[1]) * (pointEntranceXRT[1]))
      phi_radius = arctan2(-pointEntranceXRT[1],(pointEntranceXRT[0]+d)) #arccos((pointEntranceXRT[1]+d) / radius1)
      alpha = arctan(radius1 / RAYTRACER_FOCAL_LENGTH_XRT)
    vectorAfterXRTCircular[0] = radius1
    vectorAfterXRTCircular[1] = phi_radius
    vectorAfterXRTCircular[2] = alpha

    var 
      dist : seq[float]
      r1 = 0.0
      r2 = 0.0
      r3 = 0.0
      r4 = 0.0
      r5 = 0.0
      beta = 0.0 ## in degree
      r1Zyl = 0.0
      xSep = 0.0
    
    for j in 0..<allR1.len:
      if vectorAfterXRTCircular[0] > allR1[j] and vectorAfterXRTCircular[0] < allR1[j] + 0.2:
        continue
      if allR1[j] - vectorAfterXRTCircular[0] > 0.0:
        dist.add(allR1[j] - vectorAfterXRTCircular[0])
    for k in 0..<allR1.len:
      if min(dist) == allR1[k] - vectorAfterXRTCircular[0]:
        r1 = allR1[k]
        beta = degToRad(allAngles[k])
        xSep = allXsep[k]
        r2 = r1 - lMirror * sin(beta) #225mm is the length of the mirrors
        r3 = r2 - 0.5 * xSep * tan(beta)
        r4 = r3 - 0.5 * xSep * tan(3.0 * beta)
        r5 = r4 - lMirror * sin(3.0 * beta)
    if r1 == allR1[0]: continue

    var pointEntranceXRTZylKart = vec3(0.0)
    pointEntranceXRTZylKart[0] = pointEntranceXRT[0] + d 
    pointEntranceXRTZylKart[1] = pointEntranceXRT[1]
    pointEntranceXRTZylKart[2] = pointEntranceXRT[2] - centerExitPipeVT3XRT[2] 

    var pointExitCBZylKart = vec3(0.0)
    pointExitCBZylKart[0] = pointExitCB[0] + d 
    pointExitCBZylKart[1] = pointExitCB[1]
    pointExitCBZylKart[2] = pointExitCB[2] - centerExitPipeVT3XRT[2] 

    var beta3 = 3.0 * beta
    var distanceMirrors = cos(beta3) * (xSep + lMirror)
    var pointMirror1 = findPosXRT(pointEntranceXRTZylKart, pointExitCBZylKart, r1, r2, beta, lMirror, 0.0,  0.001, 1.0, 1.1)
    
    
    var s = getVectoraAfterMirror(pointEntranceXRTZylKart, pointExitCBZylKart, pointMirror1, beta, "angle") ##poinExitPipeVT3 = pointEntranceXRT
    
    
    var vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRTZylKart, pointExitCBZylKart, pointMirror1, beta, "vectorAfter")
    var pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    var pointMirror2 = findPosXRT(pointAfterMirror1, pointMirror1, r4, r5, beta3, lMirror, distanceMirrors, 0.01, 0.0, 2.5)
    var t = getVectoraAfterMirror(pointAfterMirror1, pointMirror1, pointMirror2, beta3, "angle")
    if pointMirror2[0] == 0.0 and pointMirror2[1] == 0.0 and pointMirror2[2] == 0.0: continue ## with more uncertainty, 10% of the 0.1% we loose here can be recovered, but it gets more uncertain
    
    
    #echo 2.0 * radToDeg(beta) - 0.5 * s[0]
    #echo 0.5 * t[0]
    var vectorAfterMirrors = getVectoraAfterMirror(pointAfterMirror1, pointMirror1, pointMirror2, beta3, "vectorAfter")
    

    
    #if not isDigit(intToStr(t[0].int)):
      #echo 2.0 * radToDeg(beta) - s[0]
      #echo t[0]  ## those two are the same!!! it works!!!
      #echo pointMirror1
      #echo vectorAfterMirror1
      #echo pointMirror2
      #echo vectorAfterMirrors
      #echo pointMirror1
      #ffff.add(1.1)

    ############################################# Mirrors end #################################################
    ## now get the points in the focal / detector plane

    ## because the detector is tuned in regards to the coldbore because it follows the direction of the telescope, set the origin to the detector window and turn the coordinate syste, to the detector
    var distDet = distanceMirrors - 0.5 * allXsep[8] * cos(beta) + RAYTRACER_FOCAL_LENGTH_XRT - RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW
    var epsilon = arctan(d / distDet)
    #echo epsilon
    
    var n = (distDet - pointMirror2[2]) / vectorAfterMirrors[2]

    # echo n
    var pointDetectorWindow = pointMirror2 + n * vectorAfterMirrors
    #echo pointDetectorWindow



    #for i in 1 .. < allR1.len:
      #echo lineIntersectsArea( (allR1[i] + 0.2), allR1[i], pointEntranceXRT)
      #if lineIntersectsArea( (allR1[i] + 0.2), allR1[i], pointEntranceXRT):continue

    #echo pointEntranceXRT
    var valuesPix = getPixelValue(pointEntranceXRT)
    pointdataXBefore.add(pointEntranceXRT[0])
    pointdataYBefore.add(pointEntranceXRT[1])
    pixvalsX.add(valuesPix[0])
    pixvalsY.add(valuesPix[1])

    ## Old XRT stuff ##
        #[var
          angle = (arccos(vectorBeforeXRT[2]/sqrt(vectorBeforeXRT[0]*vectorBeforeXRT[0]+vectorBeforeXRT[1]*vectorBeforeXRT[1]+vectorBeforeXRT[2]*vectorBeforeXRT[2])))  # Here we want to adress theta, the polar angle, which should be the second entrance of the vector
          r_x = pointEntranceXRT[0]
          r_y = pointEntranceXRT[1]
          theta_x = vectorBeforeXRT[0] / vectorBeforeXRT[2]
          theta_y = vectorBeforeXRT[1] / vectorBeforeXRT[2]
          theta_x_prime = theta_x - ( r_x / RAYTRACER_FOCAL_LENGTH_XRT)
          theta_y_prime = theta_y - ( r_y / RAYTRACER_FOCAL_LENGTH_XRT)]#

        #[var vectorAfterXRT = vec3(0.0)
        vectorAfterXRT[0] = vectorBeforeXRT[0] #sin(theta_x_prime) * 100.0 # theta_x_prime seemes to be in rad, since sin in c++ also does calculate from rad
        vectorAfterXRT[1] = vectorBeforeXRT[1] #sin(theta_y_prime) * 100.0
        vectorAfterXRT[2] = 100.0]#
    var centerDetectorWindow = vec3(0.0)
    centerDetectorWindow[0] = 0.0
    centerDetectorWindow[1] = 0.0
    centerDetectorWindow[2] = 0.0 #RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW # + pointEntranceXRT[2] + RAYTRACER_FOCAL_LENGTH_XRT

    var pointDetectorWindowCircle = vec3(0.0)
    pointDetectorWindowCircle[0] = tan(vectorAfterXRTCircular[2]) * RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW
    pointDetectorWindowCircle[1] = vectorAfterXRTCircular[1]
    pointDetectorWindowCircle[2] = vectorAfterXRTCircular[2]
    var pointDetectorWindow2 = vec3(0.0)
    pointDetectorWindow2[1] = pointDetectorWindowCircle[0] * sin(pointDetectorWindowCircle[1])
    pointDetectorWindow2[0] = pointDetectorWindowCircle[0] * cos(pointDetectorWindowCircle[1])
    pointDetectorWindow2[2] = 0.0


        #[var lambda_0 = ( centerDetectorWindow[2] - pointEntranceXRT[2] ) / vectorAfterXRT[2]

        pointDetectorWindow = pointEntranceXRT + lambda_0 * vectorAfterXRT
        pointDetectorWindow = pointDetectorWindow - misalignmentDetector]#

    vectorBeforeXRT = - vectorBeforeXRT # because we have to get the angles from the perspective of the XRT

    var vectorBeforeXRTPolar = vec3(0.0)  #(r,theta,phi)
    vectorBeforeXRTPolar[0] = sqrt(vectorBeforeXRT[0]*vectorBeforeXRT[0]+vectorBeforeXRT[1]*vectorBeforeXRT[1]+vectorBeforeXRT[2]*vectorBeforeXRT[2])
    vectorBeforeXRTPolar[1] = radToDeg(arccos(vectorBeforeXRT[0]/vectorBeforeXRTPolar[0]))
    vectorBeforeXRTPolar[2] = radToDeg(arctan2(vectorBeforeXRT[2],vectorBeforeXRT[1]))

    vectorBeforeXRTPolar[1] = 90.0 - vectorBeforeXRTPolar[1] #this is the pitch angle
    var p = vectorBeforeXRTPolar[1]
    vectorBeforeXRTPolar[2] = 90.0 - vectorBeforeXRTPolar[2] #this is the yaw angle, floor to roof
    var ya = vectorBeforeXRTPolar[2]
    var
      transmissionTelescopePitch = (0.0008*p*p*p*p + 1e-04*p*p*p - 0.4489*p*p - 0.3116*p + 96.787) / 100.0
      transmissionTelescopeYaw = (-0.0316*ya*ya + 0.0421*ya + 99.771) / 100.0
      probConversionMagnet = 0.025 * 9.0 * 1e-24 * (1 / (1.44 * 1.2398)) * (1 / (1.44 * 1.2398)) * (pathCB * 1e-3) * (pathCB * 1e-3)
      transmissionMagnet = cos(ya) * probConversionMagnet#1.0 # this is the transformation probability of an axion into a photon, if an axion flying straight through the magnet had one of 100%, angular dependency of the primakoff effect
      transmissionTelescopeEnergy : float
    #echo p


    if energyAx < 2.0:
      transmissionTelescopeEnergy = (-0.013086145 * energyAx * energyAx * energyAx * energyAx + 0.250552655 * energyAx * energyAx * energyAx - 1.541426299 * energyAx * energyAx + 2.064933639 * energyAx + 7.625254445) / 14.38338 #total eff area of telescope = 1438.338mm² = 14.38338cm²
    elif energyAx >= 2.0 and energyAx < 8.0:
      transmissionTelescopeEnergy = (0.0084904 * energyAx * energyAx * energyAx * energyAx * energyAx * energyAx - 0.199553 * energyAx * energyAx * energyAx * energyAx * energyAx + 1.75302 * energyAx * energyAx * energyAx * energyAx - 7.05939 * energyAx * energyAx * energyAx + 12.6706 * energyAx * energyAx - 9.23947 * energyAx + 9.96953) / 14.38338
    else:
      transmissionTelescopeEnergy = 0.0


    if transmissionTelescopeEnergy < 0.0 :
      transmissionTelescopeEnergy = 0.0

    var prob = ((43/2) * (43/2) * 3.1415) / (4 * 3.1415 * centersun[2] * centersun[2])

    var weight = ( transmissionTelescopeEnergy* transmissionTelescopePitch*transmissionTelescopeYaw* transmissionMagnet * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) ) #* prob # prob too small#transmission probabilities times time the axion spend in the magnet

    integralTotal = integralTotal + weight
    ###detector COS has (0/0) at the bottom left corner of the chip

    pointDetectorWindow[0] = pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X
    pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y
        #echo pointDetectorWindow
    pointdataX.add(pointDetectorWindow[0])
    pointdataY.add(pointDetectorWindow[1])
    weights.add(weight)


    var
      gold = ( (pointDetectorWindow[0] >= CHIPREGIONS_GOLD_X_MIN) and (pointDetectorWindow[0] <= CHIPREGIONS_GOLD_X_MAX) and (pointDetectorWindow[1] >= CHIPREGIONS_GOLD_Y_MIN) and (pointDetectorWindow[1] <= CHIPREGIONS_GOLD_Y_MAX) )
      r_xy = sqrt( ( (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X) * (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X) ) + ( (pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y) * (pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y) ) )
      silver = (r_xy <= CHIPREGIONS_SILVER_RADIUS_MAX ) and not gold
      bronze = not gold and not silver and ( r_xy <= CHIPREGIONS_BRONZE_RADIUS_MAX )
      withinWindow = r_xy < detectorWindowAperture/2
      detector = ( (pointDetectorWindow[0] >= CHIPREGIONS_CHIP_X_MIN) and (pointDetectorWindow[0] <= CHIPREGIONS_CHIP_X_MAX) and (pointDetectorWindow[1] >= CHIPREGIONS_CHIP_Y_MIN) and (pointDetectorWindow[1] <= CHIPREGIONS_CHIP_Y_MAX) )

    if(gold and withinWindow): integralGold = integralGold + weight
    if(silver and withinWindow): integralSilver = integralSilver + weight
    if(bronze and withinWindow): integralBronze = integralBronze + weight
    if(detector and withinWindow): integralDetector = integralDetector + weight


    #fluxfractionen.add(fluxfractionAll)
    #energies.add(energyAx)
    #fluxfractionen.add(1.0 + (iSun.float * 0.001))


  ## get the heatmaps out of the sequences of data X and data Y, first for the amount of data in one pixel ##
  ## compared to the overall amount and then the data in one pixel compared to the maximal amount of data in any pixel ##
  var
    beginX = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endX = 14.0 #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endY = 14.0 #- distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareheatmap(3000,3000,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,integralNormalisation)#colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareheatmap(3000,3000,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,1.0)
  #echo heatmaptable2 #= 5417.0
  echo getMaxVal(heatmaptable2, 3000)
  var heatmaptable3 = prepareheatmap(3000,3000,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,getMaxVal(heatmaptable2, 3000)) # if change number of rows: has to be in the maxVal as well
  echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
  echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]

  echo getMaxVal(heatmaptable2, 3000)
  #echo drawfancydiagrams("AxionModelFluxfraction", heatmaptable2, 3000)
  #echo drawfancydiagrams("AxionModelProbability", heatmaptable3, 3000) #Probabilities, that a photon, that hits a certain pixel could originate from an Axion, if the highest is 100%
  echo integralNormalisation # number of hits before the setup
  echo pointdataX.len # number of hits after the setup

  # get the heatmap of the data of a run for comparison
  echo "here"
  echo weights.len


  #let FILE = "likelihood_2018_2_all.h5"
  let h5file = "likelihood_2018_2.h5"
  var
    dataValuesX: seq[float]
    dataValuesY: seq[float]
  if fileExists(h5file):
    dataValuesX = getXandY(h5file,"likelihood",240, 306, "chip_3","X")
    dataValuesY = getXandY(h5file,"likelihood",240, 306, "chip_3","Y")

  var weightData : seq[float]

  var weightProb : seq[float]
  echo getLenght(heatmaptable3, 3000)
  #echo dataValuesX.len
  for i in 0 ..< dataValuesX.len:
    weightData.add(1.0)
    var X = int(floor(dataValuesX[i]*100.0))
    var Y = int(floor(dataValuesY[i]*100.0))
    #echo heatmaptable3[X][Y] * 10.0
    weightProb.add(heatmaptable3[X][Y] * 100.0)

  #echo dataValuesX
  #echo dataValuesY
  var heatmaptable4 = prepareheatmap(3000,3000,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightData,1.0)
  #echo drawfancydiagrams("AxionModelDataRaw", heatmaptable4, 3000) # the normal data of a run of chip 3
  var heatmaptable5 = prepareheatmap(3000,3000,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightProb,1.0)
  #echo drawfancydiagrams("AxionModelProbability in %", heatmaptable5, 3000) #the probability distribution of being an axion of the data of a run of chip 3
  #[var weightData1 : seq[float]
  for i in 0 ..< circleX.len:
    weightData1.add(1.0)

  var heatmaptable6 = prepareheatmap(140,140,0.0,33.0,0.0,33.0,circleX,circleY,weightData1,1.0)
  echo drawfancydiagrams("Mirrors", heatmaptable6, 140)]#

  var weightData2 : seq[float]
  for i in 0 ..< pixvalsY.len:
    weightData2.add(1.0)

  #var heatmaptable7 = prepareheatmap(1400,1400,0.0,1400.0,0.0,1400.0,pixvalsX,pixvalsY,weightData2,1.0)
  #echo drawfancydiagrams("Mirrors", heatmaptable7, 1400)
  #echo drawfancydiagrams("Mirrors",circleTotal , 1400)

  ##################################  Graphs  #######################################################

  #echo drawgraph("Energy142", radii, fluxfractionrad, "radius")
  #echo drawgraph("Radius 0%, 10% and 20%", energies, fluxfractionen, "energy") #Biased random energy distribution
  #echo drawgraph("RadiusInt", energies, fluxfracints, "energy")

  #echo drawgraph("EmissionsRates",emrates, emratesGl, "energy")
  ##echo drawgraph("EmissionsRates",energiesGl, emratesGl, "energy")
  #echo drawgraph("EmissionsRates",energies, emrates, "energy")

  #### now let's get the fluxfraction of the gold region by getting the weight of each event (probability of transition, dependend on the XRay-telescope transmission and the lenth of the
  # path the particle would have traveled through the magnet) and divide it through the number of all events


  fluxFractionTotal    = integralTotal    / integralNormalisation
  fluxFractionDetector = integralDetector / integralNormalisation
  fluxFractionBronze   = integralBronze   / integralNormalisation
  fluxFractionSilver   = integralSilver   / integralNormalisation
  fluxFractionGold     = integralGold     / integralNormalisation

  echo "Flux fraction for the gold region:"
  echo getFluxFraction("region: gold")  
  
  echo "whooop"
  echo ffff.len 





var radiationCharacteristic : string ##axionRadiation::characteristic radiationCharacteristic(axionRadiation::characteristic::sar);
radiationCharacteristic = "axionRadiation::characteristic::sar"
var coldboreBlockedLength : float64
coldboreBlockedLength = 0.0
var detectorWindowAperture : float64
detectorWindowAperture = 14.0 #mm

echo calculateFluxFractions(radiationCharacteristic, detectorWindowAperture, 0.0, 0.0, 0.0, 0.0, coldboreBlockedLength) # radiationCharacteristic = "axionRadiation::characteristic::sar"

#type
#    vector3 = ref object of Vec3

## things changed##

# weight (telescopetransmission)
# VT4 -> VT3
# XRT Focal length 1600.0 -> 1500.0
# RAYTRACER_RADIUS_PIPE_CB_VT3 = 33.6 #mm from drawing #30.0 #mm (original)
# RAYTRACER_LENGTH_PIPE_VT3_XRT = 264.7 #mm from picture #198.2 #mm (original)
# RAYTRACER_RADIUS_PIPE_VT3_XRT = 25.0 #mm from drawing #35.0 #m (original)
# 265x265 Pixel
