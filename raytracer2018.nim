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
import weave

##################rayTracer###############################

#degToRad(angle has to be done in Raytracer2014 for cos and sin

type
  ExperimentSetupKind = enum
    esCAST, esBabyIAXO

  CenterVectors = ref object
    centerEntranceCB: Vec3[float]
    centerExitCB: Vec3[float]
    centerExitPipeCBVT3: Vec3[float]
    centerExitPipeVT3XRT: Vec3[float]
    centerExitCBMagneticField: Vec3[float]
    centerSun: Vec3[float]

  ExperimentSetup* = ref object
    radiusCB*: float
    RAYTRACER_LENGTH_COLDBORE*: float
    RAYTRACER_LENGTH_COLDBORE_9T*: float
    RAYTRACER_LENGTH_PIPE_CB_VT3*: float
    radiusPipeCBVT3*: float
    RAYTRACER_LENGTH_PIPE_VT3_XRT*: float
    radiusPipeVT3XRT*: float
    RAYTRACER_FOCAL_LENGTH_XRT*: float
    distanceCBAxisXRTAxis*: float
    RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW*: float
    pipes_turned*: float
    allR1*: seq[float]
    allXsep*: seq[float]
    allAngles*: seq[float]
    lMirror*: float
    d*: float
    B*: float
    tGas*: float
    depthDet*: float

  MaterialKind = enum
    mkSi3N4 = "Si3N4"
    mkSi = "Si"
    mkAr = "Ar"

  Axion = object
    passed: bool # indicates whether axion reached the detector
    pointdataX: float
    pointdataY: float
    pointdataXBefore: float
    pointdataYBefore: float
    pointdataR: float
    weights: float
    pixvalsX: float
    pixvalsY: float
    radii: float
    energiesAx: float
    energiesAxAll: float
    kinds: MaterialKind
    transProbWindow: float
    transProbArgon: float
    transProbDetector: float
    transProbMagnet: float
    deviationDet: float
    shellNumber: int
    energiesPre: float

################################
# VARIABLES from rayTracer.h
const
  RAYTRACER_DISTANCE_SUN_EARTH = 1.5e14 #mm #ok
  radiusSun = 6.9e11                    #mm #ok
  numberOfPointsSun = 1_000_000            #100000 for statistics   

  pressGas = 14.3345 #for example P = 14.3345 mbar (corresponds to 1 bar at room temperature).


  mAxion = 0.4                          #eV for example
  g_agamma = 1e-12
var fluxesRs = newSeq[float](1000)
const
  IgnoreDetWindow = false
  IgnoreGasAbs = false

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

randomize(299792458)

func conversionProb*(B, g_agamma, length: float): float {.inline.} =
  result = 0.025 * B * B * g_agamma * g_agamma * (1 / (1.44 * 1.2398)) *
    (1 / (1.44 * 1.2398)) * (length * 1e-3) * (length * 1e-3) #g_agamma= 1e-12

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


proc getFluxFraction(chipRegionstring: string): float64 =
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
    result = getFluxFractionGold() + getFluxFractionSilver() +
        getFluxFractionBronze()

  else: echo "Error: Unknown chip region!"

#### some functions for later########
# Now let's get some functions for a random point on a disk, a biased random
# point from the solar model and a biased random energy

proc getRandomPointOnDisk(center: Vec3, radius: float64): Vec3 =
  ## This function gets a random point on a disk --> in this case this would
  ## be the exit of the coldbore ##
  var
    x = 0.0
    y = 0.0
    r = radius * sqrt(rand(1.0))
    angle = 360 * rand(1.0)
  x = cos(degToRad(angle)) * r
  y = sin(degToRad(angle)) * r
  result = vec3(x, y, 0.0) + center

proc getRandomPointFromSolarModel(center: Vec3, radius: float64,
                                  emRateCumSum: seq[float]): Vec3 =
  ## This function gives the coordinates of a random point in the sun, biased
  ## by the emissionrates (which depend on the radius and the energy) ##
  ##
  ## `emRateCumSum` is the normalized (to 1.0) cumulative sum of the emission
  ## rates over all energies at each radius
  let
    angle1 = 360 * rand(1.0)
    angle2 = 180 * rand(1.0)
    ## random number from 0 to 1 corresponding to possible solar radii.
    randEmRate = rand(1.0)
    rIdx = emRateCumSum.lowerBound(randEmRate)
    r = (0.0015 + (rIdx).float * 0.0005) * radius
  let x = cos(degToRad(angle1)) * sin(degToRad(angle2)) * r
  let y = sin(degToRad(angle1)) * sin(degToRad(angle2)) * r
  let z = cos(degToRad(angle2)) * r
  result = vec3(x, y, z) + center

proc getRandomEnergyFromSolarModel(vectorInSun, center: Vec3, radius: float64,
                                   energies: seq[float],
                                   emissionRates: seq[seq[float]],
                                   emRateCDFs: seq[seq[float]],
                                   computeType: string): float =

  ## This function gives a random energy for an event at a given radius, biased
  ## by the emissionrates at that radius. This only works if the energies to
  ## choose from are evenly distributed ##
  var
    rad = (vectorInSun - center).length
    r = rad / radius
    iRad: int
    indexRad = (r - 0.0015) / 0.0005
  if indexRad - 0.5 > floor(indexRad):
    iRad = int(ceil(indexRad))
  else: iRad = int(floor(indexRad))
  # get the normalized (to 1) CDF for this radius
  let cdfEmRate = emRateCDFs[iRad]
  # sample an index based on this CDF
  let idx = cdfEmRate.lowerBound(rand(1.0))
  ## TODO: this is horrible. Returning a physically different function from the same
  ## proc is extremely confusing! Especially given that we even calculate both at the
  ## same time!
  case computeType
  of "energy":
    let energy = energies[idx] * 0.001
    result = energy
  of "emissionRate":
    let emissionRate = emissionRates[iRad][idx]
    result = emissionRate

## The following are some functions to determine inetersection of the rays with the
## geometry (pipes, magnet, mirrors) of the setup ##
proc lineIntersectsCircle(point_1, point_2, center: Vec3,
                          radius: float64, intersect: Vec3): bool =
  ## Now a function to see if the lines from the sun will actually intersect
  ## the circle area from the magnet entrance (called coldbore) etc. ##
  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  var intersect = vec3(0.0)
  intersect = point_1 + lambda1 * vector
  let r_xy_intersect = sqrt(intersect[0] * intersect[0] + intersect[1] *
      intersect[1])
  result = r_xy_intersect < radius

proc lineIntersectsCylinderOnce(point_1: Vec3, point_2: Vec3, centerBegin: Vec3,
    centerEnd: Vec3, radius: float64, intersect: Vec3): bool =
  ## Also a function to know if the line intersected at least the whole magnet,
  ## and then only once, because else the axions would have just flown through ##
  let
    vector = point_2 - point_1
    lambda_dummy = (-1000.0 - point_1[2]) / vector[2]
    dummy = point_1 + lambda_dummy * vector
    vector_dummy = point_2 - dummy
    factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
    p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
    q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - radius*radius) / factor
    lambda_1 = -p/2.0 + sqrt(p*p/4.0 - q)
    lambda_2 = -p/2.0 - sqrt(p*p/4.0 - q)
    intersect_1 = dummy + lambda_1 * vector_dummy
    intersect_2 = dummy + lambda_2 * vector_dummy
    intersect_1_valid = (intersect_1[2] > centerBegin[2]) and
                        (intersect_1[2] < centerEnd[2])
    intersect_2_valid = (intersect_2[2] > centerBegin[2]) and
                        (intersect_2[2] < centerEnd[2])
  if ( (intersect_1_valid and intersect_2_valid) or (not intersect_1_valid and
      not intersect_2_valid)):
    return false
  elif (intersect_1_valid):
    #intersect = intersect_1
    return true
  else:
    #intersect = intersect_2
    return true

proc getIntersectLineIntersectsCylinderOnce(point_1: Vec3, point_2: Vec3,
    centerBegin: Vec3, centerEnd: Vec3, radius: float64,
    intersect: Vec3): Vec3 =
  let
    vector = point_2 - point_1
    lambda_dummy = (-1000.0 - point_1[2]) / vector[2]
    dummy = point_1 + lambda_dummy * vector
    vector_dummy = point_2 - dummy
    factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
    p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
    q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - radius*radius) / factor
    lambda_1 = -p/2.0 + sqrt(p*p/4.0 - q)
    lambda_2 = -p/2.0 - sqrt(p*p/4.0 - q)
    intersect_1 = dummy + lambda_1 * vector_dummy
    intersect_2 = dummy + lambda_2 * vector_dummy
    intersect_1_valid = (intersect_1[2] > centerBegin[2]) and
                        (intersect_1[2] < centerEnd[2])
    intersect_2_valid = (intersect_2[2] > centerBegin[2]) and
                        (intersect_2[2] < centerEnd[2])
  result = if (intersect_1_valid): intersect_1 else: intersect_2

proc circleEdges(theR1: seq[float]): seq[seq[float]] =
  ## A function to create the coordinates of the mirrors of the telescope
  ## (they're made of glass) to substract from the overall picture, because
  ## Xrays in that range don't pass through glass ##
  const
    sizeViewfield = 48.0 #mm
    d = 83.0             #mm

  var
    circleEdge = vec3(0.0)
    circleEdges = newSeqWith(1401, newSeq[float64](1401))

  for R1 in theR1:
    for phi in 34500 ..< 37500:
      for dgl in 0 .. 20:
        circleEdge[1] = (R1 + (float(dgl)*0.01)) * sin(degToRad(float(phi)*0.01))
        circleEdge[0] = (R1 + (float(dgl)*0.01)) * cos(degToRad(float(
            phi)*0.01)) - d
        var coord_X = floor(circleEdge[0] / (sizeViewfield/1400.0)) + 700
        var coord_Y = floor(circleEdge[1] / (sizeViewfield/1400.0)) + 700
        if coord_X >= 0.0:
          if coord_Y >= 0.0:
            if coord_X <= 1400.0:
              if coord_Y <= 1400.0:
                circleEdges[int(coord_X)][int(coord_Y)] = circleEdges[int(
                    coord_X)][int(coord_Y)] + 1

  result = circleEdges

proc getPixelValue(intersects: Vec3): Vec3 =
  const sizeViewfield = 48.0 #mm
  var intersectsPix = vec3(0.0)
  intersectsPix[0] = floor(intersects[0] / (sizeViewfield/1400.0)) + 700
  intersectsPix[1] = floor(intersects[1] / (sizeViewfield/1400.0)) + 700
  result = intersectsPix

proc lineIntersectsCircleEdge(circleEdges: seq[seq[float64]],
    intersectsPix: Vec3): bool =
  var
    coord_X = int(intersectsPix[0])
    coord_Y = int(intersectsPix[1])
  if circleEdges[coord_X][coord_Y] == 0.0:
    return false
  else: return true

proc lineIntersectsArea(R1: float64, prevR1: float64, intersect: Vec3): bool =
  const
    d = 83.0 #mm
  var r_xy_intersect = sqrt(intersect[1] * intersect[1] + (intersect[0]+d) * (
      intersect[0]+d))
  if r_xy_intersect <= R1 and r_xy_intersect >= (prevR1):
    return true
  else:
    return false

## Some functions to include files from outside like the run file and the emissionrate/energy files ##

proc findPosXRT*(pointXRT: Vec3, pointCB: Vec3,
                 r1, r2, angle, lMirror, distMirr, uncer, sMin, sMax: float): Vec3 =
  ## this is to find the position the ray hits the mirror shell of r1. it is after
  ## transforming the ray into a coordinate system, that has the middle of the
  ## beginning of the mirror cones as its origin
  var
    point = pointCB
    s: float
    term: float
    sMinHigh = sMin
    sMaxHigh = sMax
    pointMirror = vec3(0.0)
  let direc = pointXRT - pointCB
  template calcVal(s: float): untyped =
    let res = sqrt((point[0] + s * direc[0]) * (point[0] + s * direc[0]) + (point[
      1] + s * direc[1]) * (point[1] + s * direc[1])) - ((r2 - r1) * (point[
      2] + s * direc[2] - distMirr) / (cos(angle) * lMirror))
    res
  var mid = (sMaxHigh + sMinHigh) / 2.0
  while abs(r1 - term) > 2.0 * uncer:
    if abs(sMinHigh - sMaxHigh) < 1e-8: break
    term = calcVal(mid)
    if abs(r1 - calcVal((sMinHigh + mid) / 2.0)) < abs(r1 - calcVal((sMaxHigh + mid) / 2.0)):
      # use lower half
      sMaxHigh = mid
      mid = (sMinHigh + mid) / 2.0
    else:
      # use upper half
      sMinHigh = mid
      mid = (sMaxHigh + mid) / 2.0
  pointMirror = point + mid * direc
  result = pointMirror

proc getVectoraAfterMirror*(pointXRT, pointCB, pointMirror: Vec3,
                            angle: float, pointOrAngle: string): Vec3 =
  ## this is to find the vector after the reflection on the respective mirror
  var normalVec = vec3(0.0)
  normalVec[0] = pointMirror[0]
  normalVec[1] = pointMirror[1]
  normalVec[2] = tan(angle) * sqrt(pointMirror[0] * pointMirror[0] +
      pointMirror[1] * pointMirror[1])
  let vectorBeforeMirror = normalize(pointXRT - pointCB)
  # this is the vector product of the normal vector on pointMirror pointing
  # in the direction of the radius of the cylinder and the vector of the ray
  let vectorAxis = normalize(normalVec.cross vectorBeforeMirror)
  let alphaMirror = arcsin(abs(normalVec.dot(vectorBeforeMirror) /
                               normalVec.length))
  let vecBeforeAxis = cross(vectorBeforeMirror, vectorAxis)
  let vectorAfterMirror = vectorBeforeMirror * cos(2.0 * alphaMirror) -
      vecBeforeAxis * sin(2.0 * alphaMirror)

  let alphaTest = arccos(abs(vectorAfterMirror.dot(vectorBeforeMirror) /
      (vectorAfterMirror.length * vectorBeforeMirror.length))
  )
  var alphaVec = vec3(0.0)
  alphaVec[0] = radToDeg(alphaTest) / 2.0
  alphaVec[1] = radToDeg(alphaMirror)
  case pointOrAngle
  of "angle":
    result = alphaVec
  of "pointMirror":
    result = pointMirror
  of "pointAfter":
    result = pointMirror + 200.0 * vectorAfterMirror
  of "vectorAfter":
    result = vectorAfterMirror

proc getPointDetectorWindow(pointMirror2: Vec3, pointAfterMirror2: Vec3,
    focalLength: float, lMirror: float, xsepMiddle: float, r3Middle: float,
    dCBXray: float, pipeAngle: float): Vec3 =

  ## To calculate the point in the detector window because the pipes are turned by 3 degree (angle here in rad)
  ## First switch into new coordinate system  with its origin in the middle of the telescope and z axis turned towards the detector
  var pointMirror2Turned = vec3(0.0)
  pointMirror2Turned[0] = (pointMirror2[0] * cos(pipeAngle) + pointMirror2[2] *
      sin(pipeAngle)) - dCBXray
  pointMirror2Turned[1] = pointMirror2[1]
  pointMirror2Turned[2] = (pointMirror2[2] * cos(pipeAngle) - pointMirror2[0] *
      sin(pipeAngle)) - (lMirror + xsepMiddle / 2.0)
  var pointAfterMirror2Turned = vec3(0.0)
  pointAfterMirror2Turned[0] = (pointAfterMirror2[0] * cos(pipeAngle) +
      pointAfterMirror2[2] * sin(pipeAngle)) - dCBXray
  pointAfterMirror2Turned[1] = pointAfterMirror2[1]
  pointAfterMirror2Turned[2] = (pointAfterMirror2[2] * cos(pipeAngle) -
      pointAfterMirror2[0] * sin(pipeAngle)) - (lMirror + xsepMiddle / 2.0)
  let vectorAfterMirror2 = pointAfterMirror2Turned - pointMirror2Turned
  ## Then the distance from the middle of the telescope to the detector can be calculated with the focal length
  ## Then n can be calculated as hown many times the vector has to be applied to arrive at the detector

  var distDet = sqrt(r3Middle * r3Middle + focalLength * focalLength)
  var n = (distDet - pointMirror2Turned[2]) / vectorAfterMirror2[2]

  result = pointMirror2Turned + n * vectorAfterMirror2

proc interpTrans(fname: string): InterpolatorType[float] =
  let df = toDf(readCsv(fname, sep = ' '))
  result = newCubicSpline(df["PhotonEnergy(eV)"].toTensor(float).toRawSeq,
                          df["Transmission"].toTensor(float).toRawSeq)

# TODO: move elsewhere!
let siNfile = "./resources/Si3N4Density=3.44Thickness=0.3microns"
let spline = interpTrans(siNfile)
let siFile = "./resources/SiDensity=2.33Thickness=200.microns"
let splineStrips = interpTrans(siFile)
let detectorFile = "./resources/transmission-argon-30mm-1050mbar-295K.dat"
let splineDet = interpTrans(detectorFile)

## Now some functions for the graphs later, that store the data in heatmaps and then give them out with plotly ##

proc prepareheatmap(numberofrows: int, numberofcolumns: int,
                    start_x: float, stop_x: float, start_y: float,
                    stop_y: float,
                    data_X: seq[float], data_Y: seq[float], weight1: seq[float],
                    norm: float64): seq[seq[float]] =
  ## This function prepares a heatmap out of given X and Y values with the z value (the number of entries in a certain pixel) as the weight of the event of the X and Y value ##

  var stepsize_X = 0.0 # number of colums is the number of entries in the array in the seq and number of rows is the number of arrays in the seq
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
  # TODO: clean up!
  for i, value in data_X:
    var coord_X = floor((data_X[i] - start_x) / stepsize_X).int
    var coord_Y = floor((data_Y[i] - start_y) / stepsize_Y).int
    heatmaptable[coord_Y][coord_X] = heatmaptable[coord_Y][coord_X] + 1*weight1[i]/norm

    # if coord_X >= 0 and coord_Y >= 0 and coord_X <= float(numberofrows) and
    #   coord_Y <= float(numberofcolumns)::
    #   heatmaptable[int(coord_Y)][int(coord_X)] = heatmaptable[int(
    #     coord_Y)][int(coord_X)] + 1*weight1[i]/norm
  result = heatmaptable

proc getMaxVal(table: seq[seq[float]]): float =
  result = table.mapIt(max(it)).max

# proc serializePlot(plt: PlotJson, fname: string) =
#   ## serializes a plotly plot so that we can run it elsewhere
#   writeFile(fname, ( % plt).pretty)

proc drawfancydiagrams(diagramtitle: string,
                       objectstodraw: seq[seq[float]],
                       width: int,
                       year: string) =
  ## this function draws a hdiagram out a given heatmap ##
  var
    xs = newSeq[int](width * width)
    ys = newSeq[int](width * width)
    zs = newSeq[float](width * width)
  for y in 0 ..< width:
    for x in 0 ..< width:
      xs[y * width + x] = x
      ys[y * width + x] = y
      zs[y * width + x] = objectstodraw[y][x]

  #d.zmin = 0.0
  #d.zmax = 5e-22
  var df = seqsToDf({ "x" : xs,
                      "y" : ys,
                      "z" : zs })
    .mutate(f{float: "x-axis [mm]" ~ `x` * 14.0 / width.float},
            f{float: "y-axis [mm]" ~ `y` * 14.0 / width.float})
  template makeMinMax(knd, ax: untyped): untyped =
    template `knd ax`(): untyped =
      `CHIPREGIONS_GOLD ax knd` * width.float / 14.0
  makeMinMax(min, X)
  makeMinMax(max, X)
  makeMinMax(min, Y)
  makeMinMax(max, Y)

  echo df
  echo df.filter(f{`z` > 0.0})
  ggplot(df, aes("x-axis [mm]", "y-axis [mm]", fill = "z")) +
    geom_raster() +
    scale_x_continuous() + scale_y_continuous() + scale_fill_continuous("z") +
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = minY())) +
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = maxY())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = minX())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = maxX())) +
    ggtitle("Solar axion image for axion electron flux") +
    ggsave(&"out/axion_image_{year}.pdf")

############done with the functions, let's use them############

proc getVarsForSetup*(setup: ExperimentSetupKind): ExperimentSetup =
  # TODO: clean up, possibly make this into a toml file where one can
  # input different settings!
  result = new ExperimentSetup
  case setup
  of esCAST:
    result = ExperimentSetup(radiusCB: 21.5, #mm
      RAYTRACER_LENGTH_COLDBORE: 9756.0, #mm half B field to end of CB #ok
      RAYTRACER_LENGTH_COLDBORE_9T: 9260.0, #mm half B field to half B field #ok
      RAYTRACER_LENGTH_PIPE_CB_VT3: 2571.5, #mm should stay the same #from beam pipe drawings #ok
      radiusPipeCBVT3: 39.64, #30.0 #mm smallest aperture between end of CB and VT3
      RAYTRACER_LENGTH_PIPE_VT3_XRT: 150.0, #mm from drawings #198.2 #mm from XRT drawing #ok
      radiusPipeVT3XRT: 35.0, #25.0 #mm from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm #ok
      RAYTRACER_FOCAL_LENGTH_XRT: 1300.0, #1485.0 #mm is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 #mm was the Telescope of 2014 (MPE XRT) also: Aperatur changed #ok
      distanceCBAxisXRTAxis: 0.0, #62.1#58.44 #mm from XRT drawing #there is no difference in the axis even though the picture gets transfered 62,1mm down, but in the detector center
      RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW: 0.0, #mm #no change, because don't know
      pipes_turned: 3.0, #degree # this is the angle by which the pipes before the detector were turned in comparison to the telescope
                             # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      allR1: @[60.7095, 63.006, 65.606, 68.305, 71.105, 74.011, 77.027, 80.157,
          83.405, 86.775, 90.272, 93.902, 97.668, 101.576, 105.632], ## the radii of the shells
      allXsep: @[4.0, 4.171, 4.140, 4.221, 4.190, 4.228, 4.245, 4.288, 4.284,
          4.306, 4.324, 4.373, 4.387, 4.403, 4.481],
      #allR2: @[0.0, 60.731, 63.237, 65.838, 68.538, 71.339, 74.246, 77.263, 80.394, 83.642]
      allAngles: @[0.0, 0.579, 0.603, 0.628, 0.654, 0.680, 0.708, 0.737, 0.767,
          0.798, 0.830, 0.863, 0.898, 0.933, 0.970], ## the angles of the mirror shells coresponding to the radii above
      lMirror: 225.0, #mm Mirror length
      d: 83.0, #mm ## distance between center of colbore at XRT and center of XRT (where the focal point is on the minus x axis)
      B: 9.0, #T magnetic field of magnet
      tGas: 1.7, #K
      depthDet: 30.0 #mm
    )
  of esBabyIAXO:
    result = ExperimentSetup(radiusCB: 350.0, #mm
                             # Change:
      RAYTRACER_LENGTH_COLDBORE: 10000.0, #mm not sure if this is true
      RAYTRACER_LENGTH_COLDBORE_9T: 9990.0, #mm I know it's not 9T here should be the actual length of pipe with a stable magnetic field; can't be same length
      RAYTRACER_LENGTH_PIPE_CB_VT3: 300.0, #mm not determined
      radiusPipeCBVT3: 370.0, #mm smallest aperture between end of CB and VT4 # no Idea, I just made it wider than the coldbore
      RAYTRACER_LENGTH_PIPE_VT3_XRT: 300.0, #mm not determined
      radiusPipeVT3XRT: 370.0, #mm irrelevant, large enough to not loose anything # no idea
      RAYTRACER_FOCAL_LENGTH_XRT: 7500.0, #mm # one possibility, the other is 5050 mm
      distanceCBAxisXRTAxis: 0.0,
      RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW: 0.0, #mm #no change, because don't know #good idea
      pipes_turned: 0.0, #degree
                             # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
                             # It's again a Nustar telescope, but not just a fragment, so I suppose there are more than 13 Layers
                             # Also the widest radius in CAST was 100 mm which is not enough, it should be about 350 mm
                             # Also the angles need to be different... but I could figute them out...
      allR1: @[153.0, 154.47, 156.021, 157.652, 159.364, 161.157, 163.031, 164.985, 167.02, 169.135, 171.332, 173.609, 175.966, 178.405, 180.924, 183.524, 
          186.204, 188.965, 191.807, 194.73, 197.733, 200.817, 203.982, 207.228, 210.554, 213.961, 217.448, 221.016, 224.665, 228.395, 232.205, 236.096, 240.068, 
          244.121, 248.254, 252.468, 256.762, 261.137, 265.593, 270.13, 274.747, 279.445, 284.224, 289.084, 294.024, 299.045, 304.146, 309.329, 314.592, 319.935, 
          325.36, 330.865, 336.451, 342.117, 347.864, 353.692, 359.601, 365.59], ## the radii of the shells closest to the magnet
      allXsep: @[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      #allR2: @[0.0, 60.731, 63.237, 65.838, 68.538, 71.339, 74.246, 77.263, 80.394, 83.642]
      allAngles: @[0.29, 0.29, 0.3, 0.3, 0.3, 0.3, 0.31, 0.31, 0.32, 0.32, 0.32, 0.33, 0.33, 0.34, 0.34, 0.35, 0.35, 0.36, 0.36, 0.37, 0.37, 0.38, 0.39, 0.39, 
          0.4, 0.4, 0.41, 0.42, 0.42, 0.43, 0.44, 0.45, 0.45, 0.46, 0.47, 0.48, 0.49, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.57, 0.58, 0.59, 
          0.6, 0.61, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69], ## the angles of the mirror shells coresponding to the radii above
      lMirror: 300.0, #mm Mirror length
      d: 0.0, #mm ## distance between center of colbore at XRT and center of XRT (where the focal point is on the minus x axis)
      B: 2.0, #T magnetic field of magnet # Rather 2-3 T, not entirely homogeneous
      tGas: 293.15, #K only Gas in BabyIAXO
      depthDet: 30.0 #mm #probably not
    )

proc traceAxion(res: var Axion,
                centerVecs: CenterVectors,
                expSetup: ExperimentSetup,
                emRates: seq[seq[float]],
                emRatesRadiusCumSum: seq[float],
                emRateCDFs: seq[seq[float]],
                energies: seq[float],
                stripDistWindow, stripWidthWindow, theta: float,
                setup: ExperimentSetupKind,
                detectorWindowAperture: float,
                dfTab: Table[string, DataFrame]
               ) =
  ## Get a random point in the sun, biased by the emission rate, which is higher
  ## at smalller radii, so this will give more points in the center of the sun ##
  let pointInSun = getRandomPointFromSolarModel(centerVecs.centerSun, radiusSun, emRatesRadiusCumSum)
  #let pointInSun = getRandomPointOnDisk(centerVecs.centerSun, radiusSun)
  #let pointInSun = centerVecs.centerSun
  ## Get a random point at the end of the coldbore of the magnet to take all axions into account that make it to this point no matter where they enter the magnet ##
  let pointExitCBMagneticField = getRandomPointOnDisk(
      centerVecs.centerExitCBMagneticField, expSetup.radiusCB)

  ## Get a random energy for the axion biased by the emission rate ##
  let energyAx = getRandomEnergyFromSolarModel(
    pointInSun, centerVecs.centerSun, radiusSun, energies, emRates, emRateCDFs, "energy"
  )
  let emissionRateAx = getRandomEnergyFromSolarModel(
    pointInSun, centerVecs.centerSun, radiusSun, energies, emRates, emRateCDFs, "emissionRate"
  )
  ## Throw away all the axions, that don't make it through the piping system and therefore exit the system at some point ##
  
  # TODO: ask johanna why is `intersect` 0? Isn't being modified anywhere!
  var intersect = vec3(0.0)
  let intersectsEntranceCB = lineIntersectsCircle(pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB, expSetup.radiusCB, intersect)
  var intersectsCB = false
  res.energiesPre = energyAx
  if (not intersectsEntranceCB):
    intersectsCB = lineIntersectsCylinderOnce(pointInSun,
        pointExitCBMagneticField, centerVecs.centerEntranceCB,
        centerVecs.centerExitCBMagneticField, expSetup.radiusCB, intersect)
  if (not intersectsEntranceCB and not intersectsCB): return
  
  if (not intersectsEntranceCB): #generates problems with the weight because the weight is multiplied with the difference of the leght of the path of the particle and the legth of the coldbore
    intersect = getIntersectLineIntersectsCylinderOnce(pointInSun,
        pointExitCBMagneticField, centerVecs.centerEntranceCB,
        centerVecs.centerExitCBMagneticField, expSetup.radiusCB,
        intersect) #pointInSun + ((centerVecs.centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
  #if (not intersectsCB):
    #intersect = pointInSun + ((centerVecs.centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
  let pathCB = pointExitCBMagneticField[2] - intersect[2]
  var pointExitCB = vec3(0.0)

  if (not lineIntersectsCircle(pointInSun, pointExitCBMagneticField,
      centerVecs.centerExitCB, expSetup.radiusCB, pointExitCB)): return
  
  pointExitCB = pointInSun + ((centerVecs.centerExitCB[2] - pointInSun[2]) / (
      pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)

  var pointExitPipeCBVT3 = vec3(0.0)

  if (not lineIntersectsCircle(pointExitCBMagneticField, pointExitCB,
      centerVecs.centerExitPipeCBVT3, expSetup.radiusPipeCBVT3,
      pointExitPipeCBVT3)): 
        echo "start"
        echo "exit magnet", pointExitCBMagneticField, "exit cb", pointExitCB
        return
  
  pointExitPipeCBVT3 = pointExitCBMagneticField + ((
      centerVecs.centerExitPipeCBVT3[2] - pointExitCBMagneticField[2]) / (
      pointExitCB - pointExitCBMagneticField)[2]) * (pointExitCB - pointExitCBMagneticField)
  var pointExitPipeVT3XRT = vec3(0.0) #seq[float]

  if (not lineIntersectsCircle(pointExitCB, pointExitPipeCBVT3,
      centerVecs.centerExitPipeVT3XRT, expSetup.radiusPipeVT3XRT,
      pointExitPipeVT3XRT)): return
  
  pointExitPipeVT3XRT = pointExitCB + ((centerVecs.centerExitPipeVT3XRT[2] -
      pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) * (
      pointExitPipeCBVT3 - pointExitCB)

  var vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB
  
  ###################from the CB (coldbore(pipe in Magnet)) to the XRT (XrayTelescope)#######################

  var pointEntranceXRT = vec3(0.0)
  pointEntranceXRT[0] = pointExitPipeVT3XRT[0] #- distanceCBAxisXRTAxis
  pointEntranceXRT[1] = pointExitPipeVT3XRT[1]
  pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

  ## filters out the edges of the mirrors (NOTE: that's probably why we don't get the lines
  ## in our picture now...)
  #if (getPixelValue(pointEntranceXRT)[0] > 1400.0 or getPixelValue(pointEntranceXRT)[1] > 1400.0): return
  #if lineIntersectsCircleEdge(circleTotal, getPixelValue(pointEntranceXRT)): return

  ## there is a 2mm wide graphite block between each glass mirror, to seperate them
  ## in the middle of the X-ray telescope. Return if hit
  ## TODO: understand why this does that?!
  if pointEntranceXRT[1] <= 1.0 and pointEntranceXRT[1] >=
      -1.0: return

  ## Coordinate transform from cartesian to polar at the XRT entrance
  var
    vectorEntranceXRTCircular = vec3(0.0)
  let
    radius1 = sqrt((pointEntranceXRT[0]+expSetup.d) * (pointEntranceXRT[
        0]+expSetup.d) + (pointEntranceXRT[1]) * (pointEntranceXRT[1]))
    phi_radius = arctan2(-pointEntranceXRT[1], (pointEntranceXRT[
        0]+expSetup.d)) #arccos((pointEntranceXRT[1]+expSetup.d) / radius1)
    alpha = arctan(radius1 / expSetup.RAYTRACER_FOCAL_LENGTH_XRT)
  vectorEntranceXRTCircular[0] = radius1
  vectorEntranceXRTCircular[1] = phi_radius
  vectorEntranceXRTCircular[2] = alpha

  ## Calculate the way of the axion through the telescope by manually reflecting the ray on the two mirror layers and then ending up before the detector ##

  var
    dist: seq[float]
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    r5 = 0.0
    beta = 0.0 ## in degree
    r1Zyl = 0.0
    xSep = 0.0
    h: int

  if vectorEntranceXRTCircular[0] > expSetup.allR1[expSetup.allR1.len - 1]: return
  for j in 0..<expSetup.allR1.len:
    # get rid of where the X-rays hit the glass frontal
    if vectorEntranceXRTCircular[0] > expSetup.allR1[j] and
        vectorEntranceXRTCircular[0] < expSetup.allR1[j] + 0.2:
      return
    if expSetup.allR1[j] - vectorEntranceXRTCircular[0] > 0.0:
      dist.add(expSetup.allR1[j] - vectorEntranceXRTCircular[0])
  for k in 0..<expSetup.allR1.len:
    if min(dist) == expSetup.allR1[k] - vectorEntranceXRTCircular[0]:
      r1 = expSetup.allR1[k]
      beta = degToRad(expSetup.allAngles[k])
      xSep = expSetup.allXsep[k]
      r2 = r1 - expSetup.lMirror * sin(beta) #225mm is the length of the mirrors
      r3 = r2 - 0.5 * xSep * tan(beta)
      r4 = r3 - 0.5 * xSep * tan(3.0 * beta)
      r5 = r4 - expSetup.lMirror * sin(3.0 * beta)
      h = k
  if r1 == expSetup.allR1[0]: return

  var pointEntranceXRTZylKart = vec3(0.0)
  pointEntranceXRTZylKart[0] = pointEntranceXRT[0] + expSetup.d
  pointEntranceXRTZylKart[1] = pointEntranceXRT[1]
  pointEntranceXRTZylKart[2] = pointEntranceXRT[2] -
      centerVecs.centerExitPipeVT3XRT[2]

  var pointExitCBZylKart = vec3(0.0)
  pointExitCBZylKart[0] = pointExitCB[0] + expSetup.d
  pointExitCBZylKart[1] = pointExitCB[1]
  pointExitCBZylKart[2] = pointExitCB[2] - centerVecs.centerExitPipeVT3XRT[2]

  let
    beta3 = 3.0 * beta
    distanceMirrors = cos(beta3) * (xSep + expSetup.lMirror)
    pointMirror1 = findPosXRT(pointEntranceXRTZylKart, pointExitCBZylKart, r1,
        r2, beta, expSetup.lMirror, 0.0, 0.001, 1.0, 1.1)
  # var s = getVectoraAfterMirror(pointEntranceXRTZylKart, pointExitCBZylKart, pointMirror1, beta, "angle") ##poinExitPipeVT3 = pointEntranceXRT


  let
    vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRTZylKart,
        pointExitCBZylKart, pointMirror1, beta, "vectorAfter")
    pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    pointMirror2 = findPosXRT(pointAfterMirror1, pointMirror1, r4, r5, beta3,
        expSetup.lMirror, distanceMirrors, 0.01, 0.0, 2.5)
  # var t = getVectoraAfterMirror(pointAfterMirror1, pointMirror1, pointMirror2, beta3, "angle")
  if pointMirror2[0] == 0.0 and pointMirror2[1] == 0.0 and pointMirror2[2] ==
      0.0: return ## with more uncertainty, 10% of the 0.1% we loose here can be recovered, but it gets more uncertain
  let
    vectorAfterMirrors = getVectoraAfterMirror(pointAfterMirror1, pointMirror1,
        pointMirror2, beta3, "vectorAfter")
    pointAfterMirror2 = pointMirror2 + 200.0 * vectorAfterMirrors

  ############################################# Mirrors end #################################################
  ## now get the points in the focal / detector plane

  ## because the detector is tuned in regards to the coldbore because it follows the direction of the telescope, set the origin to the detector window and turn the coordinate syste, to the detector
  ## doesnt work for IAXO
  var
    pointDetectorWindow = vec3(0.0)
    pointEndDetector = vec3(0.0)
    pointRadialComponent: float
    distDet: float
    n: float
    n3: float
  case setup
  of esCAST:
    pointDetectorWindow = getPointDetectorWindow(pointMirror2,
        pointAfterMirror2, expSetup.RAYTRACER_FOCAL_LENGTH_XRT,
        expSetup.lMirror, expSetup.allXsep[8], 62.1, expSetup.d, degToRad(2.75))
    pointEndDetector = getPointDetectorWindow(pointMirror2, pointAfterMirror2, (
        expSetup.RAYTRACER_FOCAL_LENGTH_XRT + 30.0), expSetup.lMirror,
        expSetup.allXsep[8], 62.1, expSetup.d, degToRad(2.75))
  of esBabyIAXO:
    distDet = distanceMirrors - 0.5 * expSetup.allXsep[8] * cos(beta) +
        expSetup.RAYTRACER_FOCAL_LENGTH_XRT -
        expSetup.RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW
    n = (distDet - pointMirror2[2]) / vectorAfterMirrors[2]
    n3 = ((distDet + expSetup.depthDet) - pointMirror2[2]) / vectorAfterMirrors[2]
    pointDetectorWindow = pointMirror2 + n * vectorAfterMirrors
    pointEndDetector = pointMirror2 + n3 * vectorAfterMirrors
  #var epsilon = arctan(d / distDet)

  res.deviationDet = sqrt(pow(2.0, (pointEndDetector[0] - pointDetectorWindow[
      0])) + pow(2.0, (pointEndDetector[1] - pointDetectorWindow[1])))

  var valuesPix = getPixelValue(pointEntranceXRT)
  res.pointdataXBefore = pointEntranceXRT[0]
  res.pointdataYBefore = pointEntranceXRT[1]
  res.pixvalsX = valuesPix[0]
  res.pixvalsY = valuesPix[1]

  var centerDetectorWindow = vec3(0.0)
  centerDetectorWindow[0] = 0.0
  centerDetectorWindow[1] = 0.0
  centerDetectorWindow[2] = 0.0 #expSetup.RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW # + pointEntranceXRT[2] + expSetup.RAYTRACER_FOCAL_LENGTH_XRT

  ## Calculate the weight for each axion as the probability of its signal arriving at the detector which consists of:
  ## The probability of conversion into Xrays in the path in the homogeneous magnetic field
  ## The fraction of the flux of actions emitted from the given point in the sun per second, that actually arrives at the random point in the coldbore end
  ## The transmission probability through the LLNL Xray telescope, which depends on the energy of the Xray and its angle

  vectorBeforeXRT = - vectorBeforeXRT # because we have to get the angles from the perspective of the XRT

  var vectorBeforeXRTPolar = vec3(0.0) #(r,theta,phi)
  vectorBeforeXRTPolar[0] = vectorBeforeXRT.length
  vectorBeforeXRTPolar[1] = radToDeg(arccos(vectorBeforeXRT[
      0]/vectorBeforeXRTPolar[0]))
  vectorBeforeXRTPolar[2] = radToDeg(arctan2(vectorBeforeXRT[2], (
      vectorBeforeXRT[1])))

  #this is the pitch angle # not sure why plus 90
  vectorBeforeXRTPolar[1] = vectorBeforeXRTPolar[1] - 90.0
  let p = vectorBeforeXRTPolar[1]
  #this is the yaw angle, floor to roof
  vectorBeforeXRTPolar[2] = vectorBeforeXRTPolar[2] + 90.0
  let ya = vectorBeforeXRTPolar[2]
  let
    transmissionTelescopePitch = (0.0008*p*p*p*p + 1e-04*p*p*p - 0.4489*p*p -
        0.3116*p + 96.787) / 100.0
    # TODO: wait, is this really 6^ya etc ? Update: Also appears in thesis.
    transmissionTelescopeYaw = (6.0e-7 * pow(6.0, ya) - 1.0e-5 * pow(5.0, ya) -
        0.0001 * pow(4.0, ya) + 0.0034 * pow(3.0, ya) - 0.0292 * pow(2.0, ya) -
        0.1534 * ya + 99.959) / 100.0
    probConversionMagnet = conversionProb(expSetup.B, g_agamma, pathCB)

    distancePipe = (pointDetectorWindow[2] - pointExitCBZylKart[2]) * 1e-3 #m
    # for setup including gas
    #probConversionMagnetGas = axionConversionProb2(m_a, energyAx, pGas, tGas, (pathCB * 1e-3),
    #                                               expSetup.radiusCB, g_agamma, B)
    #absorbtionXrays = intensitySuppression2(energyAx, (pathCB * 1e-3) , distancePipe,
    #                                        pGas, tGas, roomTemp) #room temperature in K

  # this is the transformation probability of an axion into a photon, if an axion
  # flying straight through the magnet had one of 100%, angular dependency of the primakoff effect
  var transmissionMagnet: float = cos(ya) * probConversionMagnet #1.0
  # for setup with gas
  #transmissionMagnetGas = cos(ya) * probConversionMagnetGas * absorbtionXrays
  var transmissionTelescopeEnergy: float
  #echo probConversionMagnet
  case setup
  of esCAST:
    transmissionMagnet = cos(ya) * probConversionMagnet #1.0
  of esBabyIAXO:
    transmissionMagnet = cos(ya) * probConversionMagnet #0.0 #transmissionMagnetGas

  if energyAx < 2.0:
    #total eff area of telescope = 1438.338mm² = 14.38338cm² #the last thing are the mirror seperators
    transmissionTelescopeEnergy = (-0.013086145 * pow(energyAx, 4) +
                                   0.250552655 * pow(energyAx, 3) -
                                   1.541426299 * energyAx * energyAx +
                                   2.064933639 * energyAx +
                                   7.625254445) / (14.38338 - (4.3 * 0.2))
  elif energyAx >= 2.0 and energyAx < 8.0:
    transmissionTelescopeEnergy = (0.0084904 * pow(energyAx, 6) -
                                   0.199553 * pow(energyAx, 5) +
                                   1.75302 * pow(energyAx, 4) -
                                   7.05939 * pow(energyAx, 3) +
                                   12.6706 * energyAx * energyAx -
                                   9.23947 * energyAx +
                                   9.96953) / (14.38338 - (4.3 * 0.2))
  else:
    transmissionTelescopeEnergy = 0.0
  #echo "Without gas", probConversionMagnet
  #echo "With gas", probConversionMagnetGas

  if transmissionTelescopeEnergy < 0.0:
    transmissionTelescopeEnergy = 0.0
  #echo transmissionTelescopePitch

  # var prob = (expSetup.radiusCB * expSetup.radiusCB) / (4 * pointInSun[2] * pointInSun[2]) ##Flux fraction. Is this correct??* * prob * emissionRateAxPerH
  #doesnt have to be multiplied by those two as it is already time and space determined through having 1000000 axions and multiplying the total flux (produced by readOpacityFile) of axions with the entrance area of the coldbore
  # 1000000 axions that reach the coldbore then are reached after an operating time of 2.789 \times 10^{-5}\,\si{\second}
  # var emissionRateAxPerH = emissionRateAx * 3600.00
  var weight = (transmissionTelescopeEnergy *
      transmissionTelescopePitch*transmissionTelescopeYaw *
      transmissionMagnet) #transmission probabilities times axion emission rate times the flux fraction

  ##Detector window:##
  ## TODO: is the 7.0 there due to 14 mm / 2?
  if sqrt(pointDetectorWindow[0] * pointDetectorWindow[0] + pointDetectorWindow[
      1] * pointDetectorWindow[1]) > 7.0: return
  var pointDetectorWindowTurned = vec3(0.0)
  pointDetectorWindowTurned[0] = pointDetectorWindow[0] * cos(theta) +
      pointDetectorWindow[1] * sin(theta)
  pointDetectorWindowTurned[1] = pointDetectorWindow[1] * cos(theta) -
      pointDetectorWindow[0] * sin(theta)
  pointDetectorWindowTurned[2] = pointDetectorWindow[2]
  let
    x = pointDetectorWindowTurned[0]
    y = pointDetectorWindowTurned[1]

  ## Get the detector Window transmission (The stripes in the window consist of a different material than the window itself):
  var transWindow: float
  var energyAxTransWindow: int
  # TODO: assignment here of the different kinds is obviously broken. Instead of having
  # one kinds field + the others we should have some additional field or something
  ## TODO: transmission of window material etc. can also be modeled using ray tracing.
  ## probability that transmission happens at all!
  ## TODO: get the data once to avoid `toRawSeq` overhead
  if abs(y) > stripDistWindow / 2.0 and
     abs(y) < stripDistWindow / 2.0 + stripWidthWindow or
     abs(y) > 1.5 * stripDistWindow + stripWidthWindow and
     abs(y) < 1.5 * stripDistWindow + 2.0 * stripWidthWindow:
    energyAxTransWindow = dfTab["siFile"]["PhotonEnergy(eV)"].toTensor(
        float).toRawSeq.lowerBound(energyAx * 1000.0)
    transWindow = dfTab["siFile"]["Transmission"].toTensor(float)[energyAxTransWindow]
    #transWindow = splineStrips.eval(energyAx * 1000.0)
    when not IgnoreDetWindow:
      weight *= transWindow
    res.transProbWindow = transWindow
    res.transProbDetector = transWindow
    res.energiesAxAll = energyAx
    res.kinds = mkSi3N4
  else:
    energyAxTransWindow = dfTab["siNfile"]["PhotonEnergy(eV)"].toTensor(
        float).toRawSeq.lowerBound(energyAx * 1000.0)
    transWindow = dfTab["siNfile"]["Transmission"].toTensor(float)[energyAxTransWindow]
    #transWindow = spline.eval(energyAx * 1000.0)
    when not IgnoreDetWindow:
      weight *= transWindow
    res.transprobWindow = transWindow
    res.transProbDetector = transWindow
    res.energiesAxAll = energyAx
    res.kinds = mkSi

  ## Get the total probability that the Xray will be absorbed by the detector and therefore detected:
  let energyAxTransDet = dfTab["detectorFile"]["PhotonEnergy(eV)"].toTensor(
      float).toRawSeq.lowerBound(energyAx * 1000.0)
  let transDet = dfTab["detectorFile"]["Transmission"].toTensor(float)[energyAxTransWindow]
  #var transDet = splineDet.eval(energyAx * 1000.0)
  when not IgnoreGasAbs:
    weight *= 1.0 - transDet
  #echo splineDet.eval(energyAx * 1000.0)
  res.transProbArgon = transDet
  res.transProbDetector = transDet
  res.energiesAxAll = energyAx
  res.kinds = mkAr

  res.energiesAx = energyAx
  res.shellNumber = h
  ###detector COS has (0/0) at the bottom left corner of the chip
  pointRadialComponent = sqrt(pointDetectorWindow[0]*pointDetectorWindow[0]+pointDetectorWindow[1]*pointDetectorWindow[1])
  res.pointdataR = pointRadialComponent
  #echo pointRadialComponent, "x ", pointDetectorWindow[0], "y ", pointDetectorWindow[1]
  pointDetectorWindow[0] = - pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X # for the view from the detector to the sun
  pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y


  ## assign final axion position & weight
  res.pointdataX = pointDetectorWindow[0]
  res.pointdataY = pointDetectorWindow[1]
  res.weights = weight

  
  ## TODO: the following are not used for anything!
  let
    gold = ( (pointDetectorWindow[0] >= CHIPREGIONS_GOLD_X_MIN) and (
        pointDetectorWindow[0] <= CHIPREGIONS_GOLD_X_MAX) and (
        pointDetectorWindow[1] >= CHIPREGIONS_GOLD_Y_MIN) and (
        pointDetectorWindow[1] <= CHIPREGIONS_GOLD_Y_MAX))
    r_xy = sqrt( ( (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X) * (
        pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X)) + ( (
        pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y) * (
        pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y)))
    silver = (r_xy <= CHIPREGIONS_SILVER_RADIUS_MAX) and not gold
    bronze = not gold and not silver and (r_xy <= CHIPREGIONS_BRONZE_RADIUS_MAX)
    withinWindow = r_xy < detectorWindowAperture/2
    detector = ( (pointDetectorWindow[0] >= CHIPREGIONS_CHIP_X_MIN) and (
        pointDetectorWindow[0] <= CHIPREGIONS_CHIP_X_MAX) and (
        pointDetectorWindow[1] >= CHIPREGIONS_CHIP_Y_MIN) and (
        pointDetectorWindow[1] <= CHIPREGIONS_CHIP_Y_MAX))
  #echo gold
  #echo silver
  #echo bronze

  # finally set the `passed` field to indicate this axion went all the way
  res.passed = true

proc traceAxionWrapper(axBuf: ptr UncheckedArray[Axion],
                       bufLen: int,
                       centerVecs: CenterVectors,
                       expSetup: ExperimentSetup,
                       emRates: seq[seq[float]],
                       emRatesRadiusCumSum: seq[float],
                       emRateCDFs: seq[seq[float]],
                       energies: seq[float],
                       stripDistWindow, stripWidthWindow, theta: float,
                       setup: ExperimentSetupKind,
                       detectorWindowAperture: float,
                       dfTab: Table[string, DataFrame]
                      ) =
  echo "Starting weave!"
  parallelFor iSun in 0 ..< bufLen:
    captures: {axBuf, centerVecs, expSetup, emRates, emRatesRadiusCumSum, emRateCDFs,
               energies, stripDistWindow,
               stripWidthWindow, theta,
               setup, detectorWindowAperture, dfTab}
    axBuf[iSun].traceAxion(centerVecs,
                           expSetup,
                           emRates, emRatesRadiusCumSum, emRateCDFs,
                           energies,
                           stripDistWindow, stripWidthWindow, theta,
                           setup,
                           detectorWindowAperture,
                           dfTab)

proc calculateFluxFractions(axionRadiationCharacteristic: string,
                            detectorWindowAperture: float, pGas: float,
                            m_a: float, setup: ExperimentSetupKind,
                            year: string) = #The year is only for the way the window in front of the detector was turned.
  var
    distanceCBAxisXRTAxis = 0.0
    pipes_turned = 0.0
    ##Telescope transmission needs to be changed as well

  let expSetup = getVarsForSetup(setup)

  #let
  #  radiusColdbore = expSetup.radiusCB * 1e-3
  #  circleX = circleEdges( allR1)
  #  circleY = circleEdges( allR1)
  #  circleTotal = circleEdges( allR1)
  var centerSun = vec3(0.0)
  centerSun[0] = 0
  centerSun[1] = 0
  centerSun[2] = RAYTRACER_DISTANCE_SUN_EARTH

  #centerSun = centerSun + misalignmentSun

  var radiusSunTwentyPecent = radiusSun * 0.2

  var centerEntranceCB = vec3(0.0)
  centerEntranceCB[0] = 0
  centerEntranceCB[1] = 0
  centerEntranceCB[2] = 0 #coldboreBlockedLength # was 0 anyway

  var centerExitCBMagneticField = vec3(0.0)
  centerExitCBMagneticField[0] = 0
  centerExitCBMagneticField[1] = 0
  centerExitCBMagneticField[2] = expSetup.RAYTRACER_LENGTH_COLDBORE_9T

  var centerExitCB = vec3(0.0)
  centerExitCB[0] = 0
  centerExitCB[1] = 0
  centerExitCB[2] = expSetup.RAYTRACER_LENGTH_COLDBORE

  var centerExitPipeCBVT3 = vec3(0.0)
  centerExitPipeCBVT3[0] = 0
  centerExitPipeCBVT3[1] = 0
  centerExitPipeCBVT3[2] = expSetup.RAYTRACER_LENGTH_COLDBORE +
      expSetup.RAYTRACER_LENGTH_PIPE_CB_VT3

  var centerExitPipeVT3XRT = vec3(0.0)
  centerExitPipeVT3XRT[0] = 0
  centerExitPipeVT3XRT[1] = 0
  centerExitPipeVT3XRT[2] = expSetup.RAYTRACER_LENGTH_COLDBORE +
      expSetup.RAYTRACER_LENGTH_PIPE_CB_VT3 +
      expSetup.RAYTRACER_LENGTH_PIPE_VT3_XRT

  var
    integralNormalisation = 0.0
    integralTotal = 0.0
    integralDetector = 0.0
    integralBronze = 0.0
    integralSilver = 0.0
    integralGold = 0.0

  # var
    ## TODO: this doesn't belong here
    # fluxes = newSeq[float](1000)
    # energiesflux = linspace(1.0, 10000.0, 1000)

  let
    energies = linspace(1.0, 10000.0, 1112)
    roomTemp = 293.15 #K

  ## TODO: make the code use tensor for the emission rates!
  var emRatesDf = readCsv("solar_model_tensor.csv")
    .rename(f{"Radius" <- "dimension_1"}, f{"Energy" <- "dimension_2"}, f{"Flux" <- "value"})
    #.mutate(f{"Energy" ~ (`Energy` * 9 + 1.0) * 0.001})
  #ggplot(emRatesDf, aes("Radius", "Flux")) + geom_point() + ggsave("/tmp/rad_flux.pdf")
  #ggplot(emRatesDf, aes("Energy", "Flux")) + geom_point() + ggsave("/tmp/energy_flux.pdf")
  #ggplot(emRatesDf, aes("Energy", "Flux", color = factor("Radius"))) +
  #  geom_line() +
  #  xlim(0, 0.1) +
  #  ggsave("/tmp/energy_rad_flux.pdf")
  #
  #let emRatesDfRad = emRatesDf.group_by("Radius").summarize(f{float: "SumFlux" << sum(`Flux`)})
  #echo emRatesDfRad
  #ggplot(emRatesDfRad, aes("Radius", "SumFlux")) +
  #  geom_line() +
  #  ggsave("/tmp/flux_radius_sum.pdf")


  let emRatesTensor = emRatesDf["Flux"].toTensor(float)
    .reshape([emRatesDf.filter(fn {`Radius` == 0}).len, emRatesDf.filter(
        fn {`Energy` == 0}).len])
  let emRates = emRatesTensor
    .toRawSeq
    .reshape2D([emRatesTensor.shape[1], emRatesTensor.shape[0]])
  doAssert emRates[0].len == 1112
  var emRatesRadiusCumSum = emRates.mapIt(it.sum).cumSum()
  # normalize to one
  emRatesRadiusCumSum.applyIt(it / emRatesRadiusCumSum[^1])

  ## Compute all normalized CDFs of the emission rates for each radius
  var emRateCDFs = newSeq[seq[float]]()
  for iRad, f in emRates:
    var cdf = toSeq(0 ..< energies.len).mapIt(f[it] * energies[it] * energies[it] / (2 * Pi * Pi))
      .cumSum()
    cdf.applyIt(it / cdf[^1])
    emRateCDFs.add cdf

  ## sample from random point and plot
  when false:
    var es = newSeq[float]()
    var ems = newSeq[float]()
    var rs = newSeq[float]()
    var ts = newSeq[float]()
    var ps = newSeq[float]()

    for i in 0 ..< 100_000:
      let pos = getRandomPointFromSolarModel(centerSun, radiusSun, emratesRadiusCumSum)
      let r = (pos - centerSun).length()
      let energyAx = getRandomEnergyFromSolarModel(
        pos, centerSun, radiusSun, energies, emrates, emRateCDFs, "energy"
      )
      let em = getRandomEnergyFromSolarModel(
        pos, centerSun, radiusSun, energies, emrates, emRateCDFs, "emissionRate"
      )
      ts.add arccos(pos[2] / r)
      ps.add arctan(pos[1] / pos[0])
      es.add energyAx
      ems.add em
      rs.add r

    let df = seqsToDf(es, ems, rs)
    ggplot(df, aes("es")) + geom_histogram(bins = 500) + ggsave("/tmp/es.pdf")
    ggplot(df, aes("rs")) + geom_histogram(bins = 300) + ggsave("/tmp/rs.pdf")
    ggplot(df, aes("ems")) + geom_histogram(bins = 500) + ggsave("/tmp/ems.pdf")
    ggplot(df, aes("ts")) + geom_histogram() + ggsave("/tmp/ts.pdf")
    ggplot(df, aes("ps")) + geom_histogram() + ggsave("/tmp/ps.pdf")

  ## Detector Window##
  #theta angle between window strips and horizontal x axis
  var theta: float
  case year
  of "2017":
    theta = degToRad(10.8)
  of "2018":
    theta = degToRad(71.5)
  const
    stripDistWindow = 2.3  #mm
    stripWidthWindow = 0.5 #mm

  ################################################################################
  ################################################################################
  ################################################################################
  var dfTab = initTable[string, DataFrame]()
  dfTab["siFile"] = readCsv(siFile, sep = ' ')
  dfTab["siNfile"] = readCsv(siNfile, sep = ' ')
  dfTab["detectorFile"] = readCsv(detectorFile, sep = ' ')
  echo dfTab

  let centerVecs = CenterVectors(centerEntranceCB: centerEntranceCB,
                                 centerExitCB: centerExitCB,
                                 centerExitPipeCBVT3: centerExitPipeCBVT3,
                                 centerExitPipeVT3XRT: centerExitPipeVT3XRT,
                                 centerExitCBMagneticField: centerExitCBMagneticField,
                                 centerSun: centerSun)

  ## In the following we will go over a number of points in the sun, whose location and energy will be biased by the emission rate and whose track will be
  ## calculated through the CAST experimental setup from 2018 at VT3
  var axions = newSeq[Axion](numberOfPointsSun)
  var axBuf = cast[ptr UncheckedArray[Axion]](axions[0].addr)
  echo "start"
  init(Weave)
  traceAxionWrapper(axBuf, numberOfPointsSun,
                    centerVecs,
                    expSetup,
                    emRates, emRatesRadiusCumSum, emRateCDFs,
                    energies,
                    stripDistWindow, stripWidthWindow, theta,
                    setup,
                    detectorWindowAperture,
                    dfTab)
  exit(Weave)

  # walk the axions and determine `integralTotal` and `integral*`
  #if(gold and withinWindow): integralGold = integralGold + weight
  #if(silver and withinWindow): integralSilver = integralSilver + weight
  #if(bronze and withinWindow): integralBronze = integralBronze + weight
  #if(detector and withinWindow): integralDetector = integralDetector + weight

  let axionsPass = axions.filterIt(it.passed)
  echo "Passed axions ", axionsPass.len
  #echo axionsPass[0 ..< 100]
  #echo axionsPass[0 ..< 100]

  template extractPass(n: untyped): untyped =
    let n = axionsPass.mapIt(it.n)
  extractPass(deviationDet)
  extractPass(energiesAx)
  extractPass(shellNumber)

  extractPass(weights)
  extractPass(transProbArgon)
  
  extractPass(pointDataX)
  extractPass(pointDataY)
  extractPass(pointdataR)

  echo pointDataX.totensor.mean
  echo pointDataY.totensor.mean
  echo pointDataR.totensor.mean
  template extractAll(n: untyped): untyped =
    let n = axions.mapIt(it.n)
  extractAll(energiesAxAll)
  extractAll(energiesPre)
  extractAll(transProbDetector)
  extractAll(kinds)
  echo "Extracted all data!"
  ################################################################################
  ################################################################################
  ################################################################################

  let dfTransProb = seqsToDf({"Axion energy [keV]": energiesAxAll,
                               "Transmission Probability": transProbDetector,
                               "type": kinds.mapIt($it)})
  #echo energiesAxAll
  #echo transProbDetector
  #echo kinds

  # echo dfTransProb.pretty(-1)

  ggplot(dfTransProb.arrange("Axion energy [keV]"),
         aes("Axion energy [keV]", "Transmission Probability",
             color = "type")) +
    geom_line() +
    ggtitle("The transmission probability for different detector parts") +
    ggsave(&"out/TransProb_{year}.pdf")

  let dfTransProbAr = seqsToDf({"Axion energy [keV]": energiesAx,
                              "Transmission Probability": transProbArgon})
  ggplot(dfTransProbAr.arrange("Axion energy [keV]"),
         aes("Axion energy [keV]", "Transmission Probability")) +
    geom_line() +
    ggtitle("The transmission probability for the detector gas") +
    ggsave(&"out/TransProbAr_{year}.pdf")

  let dfDet = seqsToDf({"Deviation [mm]": deviationDet,
                         "Energies": energiesAx,
                         "Shell": shellNumber})
    .filter(f{Value: isNull(df["Shell"][idx]).toBool == false})

  ggplot(dfDet, aes("Deviation [mm]")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Deviation of X-rays detector entrance to readout") +
    ggsave(&"out/deviationDet_{year}.pdf")

  ggplot(dfDet, aes("Deviation [mm]", fill = factor("Shell"))) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Deviation of X-rays - detector entrance to readout") +
    ggsave(&"out/deviationDet_stacked_{year}.pdf")

  ggplot(dfDet, aes("Energies", fill = factor("Shell"))) +
    ggridges("Shell", overlap = 1.8) +
    geom_histogram(binWidth = 0.1, position = "identity") +
    ggtitle("X-ray energy distributions at detector") +
    ggsave(&"out/energies_by_shell_{year}.pdf", height = 600)

  ggplot(dfDet, aes("Deviation [mm]", fill = factor("Shell"))) +
    ggridges("Shell", overlap = 1.8) +
    geom_histogram(binWidth = 0.001, position = "identity") +
    ggtitle("Deviation of X-rays - detector entrance to readout") +
    ggsave(&"out/deviationDet_ridges_{year}.pdf", height = 600)
  
  #[let dfFluxE = seqsToDf({"Axion energy [eV]": energiesAx,
                           "Flux after experiment": fluxes})
  ggplot(dfFluxE, aes("Axion energy [eV]", "Flux after experiment")) +
    geom_point() +
    ggtitle("The flux after the experiment") +
    ggsave(&"out/FluxE_{year}.pdf")

  let dfFluxR = seqsToDf({"Axion radial position [mm]": pointdataR,
                           "Flux after experiment": fluxesR})
  ggplot(dfFluxR, aes("Axion radial position [mm]", "Flux after experiment")) +
    geom_point() +
    ggtitle("The radial flux after the experiment") +
    ggsave(&"out/FluxR_{year}.pdf")]#

  let dfRad = seqsToDf({"Radial component [mm]": pointdataR,
                        "Transmission probability": weights})
  #.filter(f{Value: isNull(df["Shell"][idx]).toBool == false})

  ggplot(dfRad, aes("Radial component [mm]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions") +
    ggsave(&"out/radialDistribution_{year}.pdf")

  let dfFluxE = seqsToDf({"Axion energy [eV]": energiesAx,
                          "Transmission probability": weights})

  ggplot(dfFluxE, aes("Axion energy [eV]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("The Axion flux after the experiment") +
    ggsave(&"out/fluxAfter_{year}.pdf") 

  let dfXY = seqsToDf({"x": pointDataX,
                       "y": pointDataY,
                       "Transmission probability": weights})
    .mutate(f{"R" ~ sqrt(`x` * `x` + `y` * `y`)})

  ggplot(dfXY, aes("x", "y")) +
    geompoint(size = some(0.5), alpha = some(0.05)) +
    ggtitle("X and Y") +
    ggsave(&"out/xy_{year}.pdf") 

  ggplot(dfXY, aes("x", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("X and Y") +
    ggsave(&"out/x_{year}.pdf")

  ggplot(dfXY, aes("R", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("R") +
    ggsave(&"out/R_{year}.pdf")

  ggplot(dfXY, aes("y", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("X and Y") +
    ggsave(&"out/y_{year}.pdf")
  
  #let dfEnergyEff = toDf(readCsv(fname2, sep = ','))
  #  .mutate(fn {"Energies [keV]" ~ `xVals` * 8.0 + 1.0},
  #          fn {"Effective Area [cm^2]" ~ `yVals` *
  #              8.5}) #.mutate(f{int -> float: "Radii" ~ `Radii`.float * 0.0005 + 0.0015})

  #ggplot(dfEnergyEff, aes("Energies [keV]", "Effective Area [cm^2]")) +
  #  geom_line() +
  #  ggtitle("The telescope energy efficiency") +
  #  ggsave(&"out/EnergyEff_{year}.pdf")

  let dfFluxE2 = seqsToDf({ "Axion energy [keV]": energiesAx,
                            "Flux after experiment": weights })
  dfFluxE2.write_csv(&"axion_gae_1e13_gagamma_{g_agamma}_flux_after_exp_N_{numberOfPointsSun}.csv")
  ggplot(dfFluxE2, aes("Axion energy [keV]", weight = "Flux after experiment")) +
    geom_histogram(binWidth = 0.1) +
    ylab("The flux after the experiment") +
    ggsave(&"out/FluxEnice_{year}.pdf")

  let dfFluxE3 = seqsToDf({"Axion energy [keV]": energiesPre})

  ggplot(dfFluxE3, aes("Axion energy [keV]")) +
    geom_histogram(binWidth = 0.1) +
    ylab("The flux before the experiment") +
    ggsave(&"out/FluxE_before_experiment_{year}.pdf")

  #[ggplot(dfFluxE, aes("Axion energy [keV]")) +#, weights = "Flux after experiment")) +
    geom_histogram() +
    ggtitle("Energy dependeny") +
    ggsave("energiesHisto.pdf")]#
  echo "all plots done, now to heatmap!"
  ## get the heatmaps out of the sequences of data X and data Y, first for the amount of data in one pixel ##
  ## compared to the overall amount and then the data in one pixel compared to the maximal amount of data in any pixel ##
  var
    beginX = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endX = 14.0  #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endY = 14.0  #- distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareheatmap(3000, 3000, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights,
      numberOfPointsSun) #colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareheatmap(256, 256, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights, 1.0)
  var heatmaptable3 = prepareheatmap(3000, 3000, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights, getMaxVal(heatmaptable2)) # if change number of rows: has to be in the maxVal as well
 # echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
 # echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]

  drawfancydiagrams("Axion Model Fluxfraction", heatmaptable2, 256, year)

  when false:
    fluxFractionTotal = integralTotal #/ integralNormalisation
    fluxFractionDetector = integralDetector #/ integralNormalisation
    fluxFractionBronze = integralBronze #/ integralNormalisation
    fluxFractionSilver = integralSilver #/ integralNormalisation
    fluxFractionGold = integralGold #/ integralNormalisation

    echo "Flux fraction for the gold region:"
    echo getFluxFraction("region: gold")
    echo "Flux fraction total"
    echo fluxFractionTotal

when isMainModule:
  # TODO: make these characteristics a mix of enums + something part of
  # `ExperimentSetup` object
  var radiationCharacteristic: string ##axionRadiation::characteristic radiationCharacteristic(axionRadiation::characteristic::sar);
  radiationCharacteristic = "axionRadiation::characteristic::sar"
  var coldboreBlockedLength: float64
  coldboreBlockedLength = 0.0
  var detectorWindowAperture: float64
  detectorWindowAperture = 14.0 #mm

  calculateFluxFractions(radiationCharacteristic, detectorWindowAperture,
                         pressGas, mAxion, esBabyIAXO,
                         "2018") # radiationCharacteristic = "axionRadiation::characteristic::sar"

  # weight (telescopetransmission)
  # VT4 -> VT3
  # XRT Focal length 1600.0 -> 1500.0
  # RAYTRACER_RADIUS_PIPE_CB_VT3 = 33.6 #mm from drawing #30.0 #mm (original)
  # RAYTRACER_LENGTH_PIPE_VT3_XRT = 264.7 #mm from picture #198.2 #mm (original)
  # RAYTRACER_RADIUS_PIPE_VT3_XRT = 25.0 #mm from drawing #35.0 #m (original)
  # 265x265 Pixel
