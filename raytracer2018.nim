# stdlib
import math, strutils, algorithm, random, sequtils, os, strformat, tables
import json except `{}`

# TODO: axionMass will be minified and the important stuff extracted
import axionMass/axionMassforMagnet

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
    allThickness*: seq[float]
    allR1*: seq[float]
    allXsep*: seq[float]
    allAngles*: seq[float]
    lMirror*: float
    d*: float
    B*: float
    pGasRoom*: float
    tGas*: float
    depthDet*: float
    radiusWindow*: float
    numberOfStrips*: int
    openAperatureRatio*: float
    windowThickness*: float
    alThickness*: float

  MaterialKind = enum
    mkSi3N4 = "Si3N4"
    mkSi = "Si"
    mkAr = "Ar"

  Axion = object
    passed: bool # indicates whether axion reached the detector
    passedTillWindow: bool
    pointdataX: float
    pointdataY: float
    pointdataXBefore: float
    pointdataYBefore: float
    pointdataR: float
    weights: float
    weightsAll: float
    transmissionMagnets: float
    yawAngles: float
    pixvalsX: float
    pixvalsY: float
    radii: float
    energiesAx: float
    energiesAxAll: float
    energiesAxWindow: float
    kinds: MaterialKind
    kindsWindow: MaterialKind
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
  numberOfPointsSun = 10_000_000            #100000 for statistics   #37734 for CAST if BabyIaxo 10 mio  #26500960 corresponding to 100_000 axions at CAST, doesnt work
  # 1000000 axions that reach the coldbore then are reached after an operating time of 2.789 \times 10^{-5}\,\si{\second} for CAST
  

  roomTemp = 293.15 #K
  mAxion = 0.0853#0.26978249412621896 #eV, corresponds to set p and T gas valus #0.4 #eV for example
  g_agamma = 1e-12

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
                          radius: float64): bool =
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

proc getIntersectlineIntersectsCircle(point_1, point_2, center: Vec3,
                          radius: float64): Vec3 =
  ## Now a function to get the intersection with one of the entrance cross sections
  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  result = point_1 + lambda1 * vector
  

proc lineIntersectsCylinderOnce(point_1: Vec3, point_2: Vec3, centerBegin: Vec3,
    centerEnd: Vec3, radius: float64): bool =
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
    centerBegin: Vec3, centerEnd: Vec3, radius: float64): Vec3 =
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


proc getPixelValue(intersects: Vec3): Vec3 =
  const sizeViewfield = 48.0 #mm
  var intersectsPix = vec3(0.0)
  intersectsPix[0] = floor(intersects[0] / (sizeViewfield/1400.0)) + 700
  intersectsPix[1] = floor(intersects[1] / (sizeViewfield/1400.0)) + 700
  result = intersectsPix


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

proc lowerBound[T](t: Tensor[T], val: T): int =
  ## returns the index of the first element in `t` that is not less than
  ## `val`, or the last element if all elements are smaller than `val`.
  ##
  ## Implementation takes the data array of the tensor, wraps it in an open
  ## array and hands it to the `algorithm.nim` procedure.
  result = lowerBound(
    toOpenArray(
      cast[ptr UncheckedArray[T]](t.unsafe_raw_offset()), 0, t.size
    ),
    val)

# proc serializePlot(plt: PlotJson, fname: string) =
#   ## serializes a plotly plot so that we can run it elsewhere
#   writeFile(fname, ( % plt).pretty)

proc drawfancydiagrams(diagramtitle: string,
                       objectstodraw: seq[seq[float]],
                       width: int,
                       year: string,
                       rSigma1: float,
                       rSigma2: float) =
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
  
  var 
    yr = linspace(- rSigma1, rSigma1, xs.len)
    yr2 = linspace(- rSigma2, rSigma2, xs.len)
  

  var df = seqsToDf({ "x" : xs,
                      "y" : ys,
                      "z" : zs,
                      "yr0": yr,
                      "yr02": yr2})
    .mutate(f{float: "x-axis [mm]" ~ `x` * 14.0 / width.float},
            f{float: "y-axis [mm]" ~ `y` * 14.0 / width.float},
            f{float: "xr" ~ sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + 7.0},
            f{float: "xrneg" ~ - sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + 7.0},
            f{float: "yr" ~ `yr0` + 7.0},
            f{float: "xr2" ~ sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + 7.0},
            f{float: "xrneg2" ~ - sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + 7.0},
            f{float: "yr2" ~ `yr02` + 7.0})
  template makeMinMax(knd, ax: untyped): untyped =
    template `knd ax`(): untyped =
      `CHIPREGIONS_GOLD ax knd` * width.float / 14.0
  makeMinMax(min, X)
  makeMinMax(max, X)
  makeMinMax(min, Y)
  makeMinMax(max, Y)

  
  
  echo df
  #echo df.filter(f{`z` > 0.0})
  ggplot(df, aes("x-axis [mm]", "y-axis [mm]", fill = "z")) +
    geom_raster() +
    scale_x_continuous() + scale_y_continuous() + scale_fill_continuous("z") +
    geompoint(aes("xr", "yr"), color = some(parseHex("F92672")), size = some(0.5)) +
    geompoint(aes("xrneg", "yr"), color = some(parseHex("F92672")), size = some(0.5)) +
    geompoint(aes("xr2", "yr2"), color = some(parseHex("ffa420")), size = some(0.5)) +
    geompoint(aes("xrneg2", "yr2"), color = some(parseHex("ffa420")), size = some(0.5)) +
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
      allThickness: @[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
          0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2], ## the radii of the shells
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
      pGasRoom: 1.0, #bar pressure of the gas
      tGas: 1.7, #K
      depthDet: 30.0, #mm
      #stripDistWindow: 2.3,  #mm
      #stripWidthWindow: 0.5, #mm #now calculated to these values
      radiusWindow: 7.0, #mm
      numberOfStrips: 4,
      openAperatureRatio: 0.838,
      windowThickness: 0.3, #microns
      alThickness: 0.02 #+-0.007
    )
  of esBabyIAXO:
    result = ExperimentSetup(radiusCB: 350.0, #mm
                             # Change:
      RAYTRACER_LENGTH_COLDBORE: 11000.0, #mm not sure if this is true but this is how its written on page 61 of the 2021 BabyIAXO paper
      RAYTRACER_LENGTH_COLDBORE_9T: 10000.0, #mm I know it's not 9T here should be the actual length of pipe with a stable magnetic field; can't be same length
      RAYTRACER_LENGTH_PIPE_CB_VT3: 300.0, #mm not determined
      radiusPipeCBVT3: 370.0, #mm smallest aperture between end of CB and VT4 # no Idea, I just made it wider than the coldbore
      RAYTRACER_LENGTH_PIPE_VT3_XRT: 300.0, #mm not determined
      radiusPipeVT3XRT: 370.0, #mm irrelevant, large enough to not loose anything # no idea
      RAYTRACER_FOCAL_LENGTH_XRT: 7500.0, #mm # one possibility, the other is 5050 mm
      distanceCBAxisXRTAxis: 0.0,
      RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW: 0.0, #mm #no change, because don't know #good idea
      pipes_turned: 0.0, #degree
                             # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      #allR3: @[151.61, 153.88, 156.17, 158.48, 160.82, 163.18, 165.57, 167.98, 170.42, 172.88, 175.37, 177.88, 180.42, 183.14, 185.89, 188.67, 191.48, 
          #194.32, 197.19, 200.09, 203.02, 206.03, 209.07, 212.14, 215.24, 218.37, 221.54, 224.74, 227.97, 231.24, 234.54, 237.87, 241.24, 244.85, 
          #248.5, 252.19, 255.92, 259.68, 263.48, 267.32, 271.2, 275.12, 279.08, 283.09, 287.14, 291.38, 295.72, 300.11, 304.54, 309.02, 313.54, 
          #318.11, 322.73, 327.4, 332.12, 336.88, 341.69, 346.55] #these are the real values but R3
      allThickness: @[0.468, 0.475, 0.482, 0.490, 0.497, 0.504, 0.511, 0.519, 0.526, 0.534, 0.542, 0.549, 0.557, 0.566, 0.574, 0.583, 0.591, 0.600, 0.609, 0.618, 
          0.627, 0.636, 0.646, 0.655, 0.665, 0.675, 0.684, 0.694, 0.704, 0.714, 0.724, 0.735, 0.745, 0.756, 0.768, 0.779, 0.790, 0.802, 0.814, 0.826, 0.838, 0.850, 
          0.862, 0.874, 0.887, 0.900, 0.913, 0.927, 0.941, 0.955, 0.968, 0.983, 0.997, 1.011, 1.026, 1.041, 1.055, 1.070],
      allR1: @[153.126, 155.419, 157.731, 160.065, 162.428, 164.812, 167.225, 169.66, 172.124, 174.608, 177.123, 179.658, 182.224, 184.971, 187.749, 190.556, 
          193.394, 196.263, 199.161, 202.09, 205.05, 208.09, 211.16, 214.261, 217.392, 220.553, 223.755, 226.987, 230.249, 233.552, 236.885, 240.248, 243.652, 
          247.298, 250.984, 254.711, 258.478, 262.276, 266.114, 269.992, 273.911, 277.87, 281.869, 285.92, 290.01, 294.292, 298.676, 303.109, 307.584, 312.108, 
          316.674, 321.289, 325.955, 330.672, 335.439, 340.246, 345.104, 350.013], ## the radii of the shells closest to the magnet, now correct
      allXsep: @[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      allAngles: @[0.29, 0.294, 0.298, 0.303, 0.307, 0.312, 0.316, 0.321, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 0.36, 0.366, 0.371, 0.377, 0.382, 0.388, 
          0.393, 0.399, 0.405, 0.411, 0.417, 0.423, 0.429, 0.435, 0.441, 0.448, 0.454, 0.461, 0.467, 0.474, 0.481, 0.489, 0.496, 0.503, 0.51, 0.518, 0.525, 0.533, 
          0.54, 0.548, 0.556, 0.564, 0.573, 0.581, 0.59, 0.598, 0.607, 0.616, 0.625, 0.634, 0.643, 0.652, 0.661], ## the angles of the mirror shells coresponding to the radii above, now correct
      lMirror: 300.0, #mm Mirror length
      d: 0.0, #mm ## distance between center of colbore at XRT and center of XRT (where the focal point is on the minus x axis)
      B: 2.0, #T magnetic field of magnet # Rather 2-3 T, not entirely homogeneous
      pGasRoom: 1.0, #bar, pressure of the gas #for example P = 14.3345 mbar (corresponds to 1 bar at room temperature).
      tGas: 100.0, #293.15, #K only Gas in BabyIAXO
      depthDet: 30.0, #mm #probably not
      radiusWindow: 4.0, #mm
      numberOfStrips: 20, #maybe baby
      openAperatureRatio: 0.95,
      windowThickness: 0.1, #microns #options are: 0.3, 0.15 and 0.1
      alThickness: 0.015
    )

proc calcWindowVals(radiusWindow: float, numberOfStrips: int, openAperatureRatio: float): seq[float] =
  let 
    totalArea = radiusWindow * radiusWindow * PI
    areaOfStrips = totalArea * (1.0 - openAperatureRatio)
    dAndwPerStrip = radiusWindow * 2.0 / (numberOfStrips.float + 1.0)  #width and distance between strips per strip; 
                                                                       #the width on both sides is a whole width 
                                                                       #(Don't know what to do about the additional string width) 
                                                                       #not important at high strip number but at low like CAST

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
    stripsWidthandDist = @[widthStrips, distStrips]
  echo widthStrips
  echo distStrips
  result = stripsWidthandDist

proc traceAxion(res: var Axion,
                centerVecs: CenterVectors,
                expSetup: ExperimentSetup,
                emRates: seq[seq[float]],
                emRatesRadiusCumSum: seq[float],
                emRateCDFs: seq[seq[float]],
                energies: seq[float],
                stripDistWindow, stripWidthWindow, theta: float,
                setup: ExperimentSetupKind,
                year: string, 
                stage: string,
                detectorWindowAperture: float,
                dfTab: Table[string, DataFrame],
                dfTable: Table[string, DataFrame]
               ) =
  ## Get a random point in the sun, biased by the emission rate, which is higher
  ## at smalller radii, so this will give more points in the center of the sun ##
  let pointInSun = getRandomPointFromSolarModel(centerVecs.centerSun, radiusSun, emRatesRadiusCumSum)

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
  
  let intersectsEntranceCB = lineIntersectsCircle(pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB, expSetup.radiusCB)
  var intersectsCB = false
  res.energiesPre = energyAx
  if (not intersectsEntranceCB):
    intersectsCB = lineIntersectsCylinderOnce(pointInSun,
        pointExitCBMagneticField, centerVecs.centerEntranceCB,
        centerVecs.centerExitCBMagneticField, expSetup.radiusCB)
  if (not intersectsEntranceCB and not intersectsCB): return

  var intersect = vec3(0.0) #isnt't changed for axions that hit the entrance of the coldbore because the z value is 0 then anyways
  if (not intersectsEntranceCB): #generates problems with the weight because the weight is multiplied with the difference of the leght of the path of the particle and the legth of the coldbore
    intersect = getIntersectLineIntersectsCylinderOnce(pointInSun,
        pointExitCBMagneticField, centerVecs.centerEntranceCB,
        centerVecs.centerExitCBMagneticField, expSetup.radiusCB) #pointInSun + ((centerVecs.centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
  else: intersect = getIntersectlineIntersectsCircle(pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB, expSetup.radiusCB)

  ##get the length of the path of the axion in the magnetic field to get the probability of conversion later
  let pathCB = (pointExitCBMagneticField - intersect).length 
  var pointExitCB = vec3(0.0)
  if (not lineIntersectsCircle(pointInSun, pointExitCBMagneticField,
      centerVecs.centerExitCB, expSetup.radiusCB)): return
  
  pointExitCB = pointInSun + ((centerVecs.centerExitCB[2] - pointInSun[2]) / (
      pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)

  var pointExitPipeCBVT3 = vec3(0.0)

  if (not lineIntersectsCircle(pointExitCBMagneticField, pointExitCB,
      centerVecs.centerExitPipeCBVT3, expSetup.radiusPipeCBVT3)): 
        echo "start"
        echo "exit magnet", pointExitCBMagneticField, "exit cb", pointExitCB
        return
  
  pointExitPipeCBVT3 = pointExitCBMagneticField + ((
      centerVecs.centerExitPipeCBVT3[2] - pointExitCBMagneticField[2]) / (
      pointExitCB - pointExitCBMagneticField)[2]) * (pointExitCB - pointExitCBMagneticField)
  var pointExitPipeVT3XRT = vec3(0.0) 

  if (not lineIntersectsCircle(pointExitCB, pointExitPipeCBVT3,
      centerVecs.centerExitPipeVT3XRT, expSetup.radiusPipeVT3XRT)): return
  
  pointExitPipeVT3XRT = pointExitCB + ((centerVecs.centerExitPipeVT3XRT[2] -
      pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) * (
      pointExitPipeCBVT3 - pointExitCB)

  var vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB
  
  ###################from the CB (coldbore(pipe in Magnet)) to the XRT (XrayTelescope)#######################

  var pointEntranceXRT = vec3(0.0)
  pointEntranceXRT[0] = pointExitPipeVT3XRT[0] #- distanceCBAxisXRTAxis
  pointEntranceXRT[1] = pointExitPipeVT3XRT[1]
  pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

  ## Coordinate transform from cartesian to polar at the XRT entrance
  var
    vectorEntranceXRTCircular = vec3(0.0)
  let
    radius1 = sqrt((pointEntranceXRT[0]+expSetup.d) * (pointEntranceXRT[
        0]+expSetup.d) + (pointEntranceXRT[1]) * (pointEntranceXRT[1]))
    phi_radius = arctan2(-pointEntranceXRT[1], (pointEntranceXRT[
        0]+expSetup.d)) #arccos((pointEntranceXRT[1]+expSetup.d) / radius1)
    alpha = arctan(radius1 / expSetup.RAYTRACER_FOCAL_LENGTH_XRT)
    phi_flat = radtoDeg(arccos(pointEntranceXRT[0] / radius1))
  vectorEntranceXRTCircular[0] = radius1
  vectorEntranceXRTCircular[1] = phi_radius #in rad
  vectorEntranceXRTCircular[2] = alpha #in rad

  
  ## there is a 2mm wide graphite block between each glass mirror, to seperate them
  ## in the middle of the X-ray telescope. Return if hit
  ## BabyIAXO return if X-ray its the spider structure
  case setup
  of esCAST:
    if pointEntranceXRT[1] <= 1.0 and pointEntranceXRT[1] >=
        -1.0: return
  of esBabyIAXO:
    ## here we have a spider structure for the XMM telescope:
    if vectorEntranceXRTCircular[0] <= 64.7: #doesnt really matter because there are no mirrors in the middle and these axions dont reach the window anyways
      return
    elif vectorEntranceXRTCircular[0] < 151.6 and vectorEntranceXRTCircular[0] > (151.6 - 20.9): #doesnt really matter because there are no mirrors in the middle and these axions dont reach the window anyways
      return
    for i in 0..16:
      if (phi_flat >= (-1.25 + 22.5 * i.float) and phi_flat <= (1.25 + 22.5 * i.float)): #spider strips (actually wider for innermost but doesnn't matter because it doesnt reach the window anyways)
        return
    #TODO: inner spider structure that doesnt matter

  
  ## Calculate the way of the axion through the telescope by manually reflecting the ray on the two mirror layers and then ending up before the detector ##

  var
    dist: seq[float]
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    r5 = 0.0
    beta = 0.0 ## in degree
    xSep = 0.0
    h: int

  if vectorEntranceXRTCircular[0] > expSetup.allR1[expSetup.allR1.len - 1]: return
  for j in 0..<expSetup.allR1.len:
    # get rid of where the X-rays hit the glass frontal
    if vectorEntranceXRTCircular[0] > expSetup.allR1[j] and
        vectorEntranceXRTCircular[0] < expSetup.allR1[j] + expSetup.allThickness[j]: 
      return
    if expSetup.allR1[j] - vectorEntranceXRTCircular[0] > 0.0:
      dist.add(expSetup.allR1[j] - vectorEntranceXRTCircular[0])
  for k in 0..<expSetup.allR1.len:
    if min(dist) == expSetup.allR1[k] - vectorEntranceXRTCircular[0]:
      r1 = expSetup.allR1[k]
      beta = degToRad(expSetup.allAngles[k])
      xSep = expSetup.allXsep[k]
      r2 = r1 - expSetup.lMirror * sin(beta) 
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


  let
    vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRTZylKart,
        pointExitCBZylKart, pointMirror1, beta, "vectorAfter")
    pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    pointMirror2 = findPosXRT(pointAfterMirror1, pointMirror1, r4, r5, beta3,
        expSetup.lMirror, distanceMirrors, 0.01, 0.0, 2.5)
  if pointMirror2[0] == 0.0 and pointMirror2[1] == 0.0 and pointMirror2[2] ==
      0.0: return ## with more uncertainty, 10% of the 0.1% we loose here can be recovered, but it gets more uncertain
  let
    vectorAfterMirrors = getVectoraAfterMirror(pointAfterMirror1, pointMirror1,
        pointMirror2, beta3, "vectorAfter")
    pointAfterMirror2 = pointMirror2 + 200.0 * vectorAfterMirrors

  ############################################# Mirrors end #################################################
  ## now get the points in the focal / detector plane

  ## because the detector is tuned in regards to the coldbore because it follows the direction of the telescope, set the origin to the detector window and 
  ## turn the coordinate syste, to the detector for CAST

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
  ## NOT! The fraction of the flux of actions emitted from the given point in the sun per second, that actually arrives at the random point in the coldbore end, because this is handled by biasing the random origin
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
    # TODO: change this for CAST in the same way BabyIAXO was done with the direct reflectivity
    transmissionTelescopeYaw = (6.0e-7 * pow(ya, 6.0) - 1.0e-5 * pow(ya, 5.0) -
        0.0001 * pow(ya, 4.0) + 0.0034 * pow(ya, 3.0) - 0.0292 * pow(ya, 2.0) -
        0.1534 * ya + 99.959) / 100.0
    probConversionMagnet = conversionProb(expSetup.B, g_agamma, pathCB)
  
    distancePipe = (pointDetectorWindow[2] - pointExitCBZylKart[2]) * 1e-3 #m
    pGas = expSetup.pGasRoom / roomTemp * expSetup.tGas
    effPhotonMass = effPhotonMass2(pGas, (pathCB * 1e-3), (expSetup.radiusCB * 1e-3), expSetup.tGas)
    probConversionMagnetGas = axionConversionProb2(mAxion, energyAx, pGas, expSetup.tGas, (pathCB * 1e-3),
                                                   (expSetup.radiusCB * 1e-3), g_agamma, expSetup.B) # for setup including gas: functions are in axionmass/axionMassforMagnet
    absorbtionXrays = intensitySuppression2(energyAx, (pathCB * 1e-3), distancePipe,
                                            pGas, expSetup.tGas, roomTemp) #room temperature in K
  #echo "axion mass in [eV] ", mAxion, " effective photon mass in [eV] ", effPhotonMass

  ## this is the transformation probability of an axion into a photon, if an axion
  ## flying straight through the magnet had one of 100%, angular dependency of the primakoff effect
  var 
    weight:float
    transmissionMagnet: float
    transmissionMagnetGas = cos(ya) * probConversionMagnetGas * absorbtionXrays # for setup with gas
    transmissionTelescopeEnergy: float
  
  
  case stage
  of "vacuum":
    transmissionMagnet = cos(ya) * probConversionMagnet 
  of "gas":
    transmissionMagnet = transmissionMagnetGas
  
  res.transmissionMagnets = transmissionMagnet
  res.yawAngles = ya

  let
    angle1 = getVectoraAfterMirror(pointEntranceXRTZylKart,
    pointExitCBZylKart, pointMirror1, beta, "angle")
    angle2 = getVectoraAfterMirror(pointAfterMirror1, pointMirror1,
        pointMirror2, beta3, "angle")
    alpha1 = angle1[1].round(2)
    alpha2 = angle2[1].round(2)

  case setup
  of esCAST:
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
    if transmissionTelescopeEnergy < 0.0:
      transmissionTelescopeEnergy = 0.0
    weight = (transmissionTelescopeEnergy *
      transmissionTelescopePitch*transmissionTelescopeYaw *
      transmissionMagnet) #transmission probabilities times axion emission rate times the flux fraction
  of esBabyIAXO:
    let
      energyAxReflection1 = dfTable[fmt"goldfile{alpha1:4.2f}"]["PhotonEnergy(eV)"].toTensor(
              float).lowerBound(energyAx * 1000.0)
      energyAxReflection2 = dfTable[fmt"goldfile{alpha2:4.2f}"]["PhotonEnergy(eV)"].toTensor(
              float).lowerBound(energyAx * 1000.0)
      reflectionProb1 = dfTable[fmt"goldfile{alpha1:4.2f}"]["Reflectivity"].toTensor(float)[energyAxReflection1]
      reflectionProb2 = dfTable[fmt"goldfile{alpha2:4.2f}"]["Reflectivity"].toTensor(float)[energyAxReflection2]
    weight = reflectionProb1 * reflectionProb2 * transmissionMagnet#also yaw and pitch dependend
  #if alpha1 < 0.45 or alpha1 > 1.11:
    #echo alpha1, " ", alpha2


  if weight != 0:
    res.passedTillWindow = true

  ##Detector window:##
  if sqrt(pointDetectorWindow[0] * pointDetectorWindow[0] + pointDetectorWindow[
      1] * pointDetectorWindow[1]) > expSetup.radiusWindow: return
  var pointDetectorWindowTurned = vec3(0.0)
  pointDetectorWindowTurned[0] = pointDetectorWindow[0] * cos(theta) +
      pointDetectorWindow[1] * sin(theta)
  pointDetectorWindowTurned[1] = pointDetectorWindow[1] * cos(theta) -
      pointDetectorWindow[0] * sin(theta)
  pointDetectorWindowTurned[2] = pointDetectorWindow[2]
  let
    x = pointDetectorWindowTurned[0]
    y = pointDetectorWindowTurned[1]

  ## Get the detector Window transmission (The stripes in the window consist of a different
  ## material than the window itself)
  var transWindow: float
  var energyAxTransWindow: int
  # TODO: assignment here of the different kinds is obviously broken. Instead of having
  # one kinds field + the others we should have some additional field or something #made two assignments and now it works
  ## TODO: transmission of window material etc. can also be modeled using ray tracing.
  ## probability that transmission happens at all!
  for i in 0..(expSetup.numberOfStrips/2).round.int - 1:
    if abs(y) > (1.0 * i.float + 0.5) * stripDistWindow + i.float * stripWidthWindow and
      abs(y) < (1.0 * i.float + 0.5) * stripDistWindow + (i.float + 1.0) * stripWidthWindow:
      energyAxTransWindow = dfTab["siFile"]["PhotonEnergy(eV)"].toTensor(
          float).lowerBound(energyAx * 1000.0)
      transWindow = dfTab["siFile"]["Transmission"].toTensor(float)[energyAxTransWindow] *
                    dfTab["alFile"]["Transmission"].toTensor(float)[energyAxTransWindow]
      when not IgnoreDetWindow:
        weight *= transWindow
      res.transProbWindow = transWindow
      res.transProbDetector = transWindow
      res.energiesAxAll = energyAx
      res.energiesAxWindow = energyAx
      res.kinds = mkSi
      res.kindsWindow = mkSi
    else:
      energyAxTransWindow = dfTab["siNfile"]["PhotonEnergy(eV)"].toTensor(
          float).lowerBound(energyAx * 1000.0)
      transWindow = dfTab["siNfile"]["Transmission"].toTensor(float)[energyAxTransWindow] *
                    dfTab["alFile"]["Transmission"].toTensor(float)[energyAxTransWindow]
      when not IgnoreDetWindow:
        weight *= transWindow
      res.transprobWindow = transWindow
      res.transProbDetector = transWindow
      res.energiesAxAll = energyAx
      res.energiesAxWindow = energyAx
      res.kinds = mkSi3N4
      res.kindsWindow = mkSi3N4

  ## Get the total probability that the Xray will be absorbed by the detector and therefore detected:
  let energyAxTransDet = dfTab["detectorFile"]["PhotonEnergy(eV)"].toTensor(
      float).lowerBound(energyAx * 1000.0)
  let transDet = dfTab["detectorFile"]["Transmission"].toTensor(float)[energyAxTransDet]

  when not IgnoreGasAbs:
    weight *= 1.0 - transDet
  res.transProbArgon = transDet
  res.transProbDetector = transDet

  res.energiesAxAll = energyAx
  res.kinds = mkAr
  res.energiesAx = energyAx
  res.shellNumber = h
  
  ###detector COS has (0/0) at the bottom left corner of the chip
  pointRadialComponent = sqrt(pointDetectorWindow[0]*pointDetectorWindow[0]+pointDetectorWindow[1]*pointDetectorWindow[1])
  res.pointdataR = pointRadialComponent
  pointDetectorWindow[0] = - pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X # for the view from the detector to the sun
  pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y


  ## assign final axion position & weight
  res.pointdataX = pointDetectorWindow[0]
  res.pointdataY = pointDetectorWindow[1]
  res.weights = weight
  res.weightsAll = weight

  
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

  ## finally set the `passed` field to indicate this axion went all the way 
  if weight != 0:
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
                       year: string, 
                       stage: string,
                       detectorWindowAperture: float,
                       dfTab: Table[string, DataFrame],
                       dfTable: Table[string, DataFrame]
                       ) =
  echo "Starting weave!"
  parallelFor iSun in 0 ..< bufLen:
    captures: {axBuf, centerVecs, expSetup, emRates, emRatesRadiusCumSum, emRateCDFs,
               energies, stripDistWindow,
               stripWidthWindow, theta,
               setup, year, stage, detectorWindowAperture, dfTab, dfTable}
    axBuf[iSun].traceAxion(centerVecs,
                           expSetup,
                           emRates, emRatesRadiusCumSum, emRateCDFs,
                           energies,
                           stripDistWindow, stripWidthWindow, theta,
                           setup,
                           year,
                           stage,
                           detectorWindowAperture,
                           dfTab, dfTable)

proc calculateFluxFractions(axionRadiationCharacteristic: string,
                            detectorWindowAperture: float,
                            setup: ExperimentSetupKind,
                            year: string,
                            stage: string) = #The year is only for the way the window in front of the detector was turned.
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

  let
    energies = linspace(1.0, 10000.0, 1112)



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

  
  
  ################################################################################
  #############################Detector Window####################################
  ################################################################################  
  var
    stripWidthWindow = calcWindowVals(expSetup.radiusWindow, expSetup.numberOfStrips, expSetup.openAperatureRatio)[0]
    stripDistWindow = calcWindowVals(expSetup.radiusWindow, expSetup.numberOfStrips, expSetup.openAperatureRatio)[1]
    theta: float   #theta angle between window strips and horizontal x axis
  case year
  of "2017":
    theta = degToRad(10.8)
  of "2018":
    theta = degToRad(71.5)

  let siNfile = &"./resources/Si3N4Density=3.44Thickness={expSetup.windowThickness}microns"
  let siFile = "./resources/SiDensity=2.33Thickness=200.microns"
  let detectorFile = "./resources/transmission-argon-30mm-1050mbar-295K.dat"
  let alFile = &"./resources/AlDensity=2.7Thickness={expSetup.alThickness}microns"


  var dfTab = initTable[string, DataFrame]()
  dfTab["siFile"] = readCsv(siFile, sep = ' ')
  dfTab["siNfile"] = readCsv(siNfile, sep = ' ')
  dfTab["detectorFile"] = readCsv(detectorFile, sep = ' ')
  dfTab["alFile"] = readCsv(alFile, sep = ' ')
  

  var 
    goldfile: string
    alpha: float
    dfTable = initTable[string, DataFrame]()

  for i in 13..83:
    alpha = (i.float * 0.01).round(2)
    goldfile = fmt"./resources/reflectivity/{alpha:4.2f}degGold0.25microns"
    dfTable[fmt"goldfile{alpha:4.2f}"] = readCsv(goldfile, sep = ' ')

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
                    year,
                    stage,
                    detectorWindowAperture,
                    dfTab,
                    dfTable)
  exit(Weave)

  # walk the axions and determine `integralTotal` and `integral*`
  #if(gold and withinWindow): integralGold = integralGold + weight
  #if(silver and withinWindow): integralSilver = integralSilver + weight
  #if(bronze and withinWindow): integralBronze = integralBronze + weight
  #if(detector and withinWindow): integralDetector = integralDetector + weight

  let axionsPass = axions.filterIt(it.passed)
  echo "Passed axions ", axionsPass.len
  let axionsPassW = axions.filterIt(it.passedTillWindow)
  echo "Passed axions until the Window ", axionsPassW.len

  template extractPass(n: untyped): untyped =
    let n = axionsPass.mapIt(it.n)  
  template extractAll(n: untyped): untyped =
    let n = axions.mapIt(it.n)
  extractPass(deviationDet)
  extractPass(energiesAx)
  extractPass(shellNumber)

  extractPass(weights)
  extractPass(transmissionMagnets)
  extractPass(yawAngles)
  extractPass(transProbArgon)

  extractPass(pointDataX)
  extractPass(pointDataY)
  extractPass(pointdataR)

  echo pointDataX.totensor.mean
  echo pointDataY.totensor.mean
  echo pointDataR.totensor.mean
  
  extractAll(weightsAll)
  extractAll(energiesAxAll)
  extractAll(energiesAxWindow)
  extractAll(energiesPre)
  extractAll(transProbDetector)
  extractAll(transprobWindow)
  extractAll(kinds)
  extractAll(kindsWindow)
  echo "Extracted all data!"
  ################################################################################
  ################################################################################
  ################################################################################

  let dfTransProb = seqsToDf({"Axion energy [keV]": energiesAxAll,
                               "Transmission Probability": transProbDetector,
                               "type": kinds.mapIt($it),
                               "Axion energy window[keV]":energiesAxWindow,
                               "Transmission Probability window": transprobWindow,
                               "type window":kindsWindow.mapIt($it),
                               "Flux after experiment": weightsAll})
  
  
  ggplot(dfTransProb.arrange("Axion energy [keV]")
         ) +
    geom_line(aes("Axion energy [keV]", "Transmission Probability",
             color = "type")) +
    geom_line(aes("Axion energy [keV]", "Transmission Probability window",
             color = "type window")) +
    geom_histogram(aes("Axion energy [keV]", weight = "Flux after experiment"), binWidth = 0.01) +
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
  

  let dfRad = seqsToDf({"Radial component [mm]": pointdataR,
                        "Transmission probability": weights,
                        "x": pointDataX,
                       "y": pointDataY})
  
  echo dfRad.arrange("Radial component [mm]")
  ggplot(dfRad, aes("Radial component [mm]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions") +
    ggsave(&"out/radialDistribution_{year}.pdf")
  
  let dfFluxE = seqsToDf({"Axion energy [keV]": energiesAx,
                          "Transmission probability": weights})

  ggplot(dfFluxE, aes("Axion energy [keV]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.0001) +
    ggtitle("The Axion flux after the experiment") +
    ggsave(&"out/fluxAfter_{year}.pdf") 


  let dfXY = seqsToDf({"x": pointDataX,
                       "y": pointDataY,
                       "Transmission probability": weights})
    .mutate(f{"R" ~ sqrt(`x` * `x` + `y` * `y`)})


  ggplot(dfXY, aes("x", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("X and Y") +
    ggsave(&"out/x_{year}.pdf")

  ggplot(dfXY, aes("R", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("R") +
    ggsave(&"out/R_{year}.pdf")

  let dfMag = seqsToDf({"Transmission probability": transmissionMagnets,
                       "Angles between path and magnetic field": yawAngles,
                       "Axion energy[keV]": energiesAx}).arrange("Angles between path and magnetic field")


  ggplot(dfMag, aes("Angles between path and magnetic field", "Transmission probability")) +
    geom_point(size = some(0.5), alpha = some(0.1)) +
    ylim(3.1e-24, 3.16e-24) +
    ggtitle("The probability of the transformation of axions to X-rays in the magnet") +
    ggsave(&"out/transMagnet_{year}.pdf")

  ggplot(dfMag, aes("Axion energy[keV]", "Transmission probability")) +
    geom_point(size = some(0.5), alpha = some(0.1)) +
    #ylim(3.1e-24, 3.16e-24) +
    ggtitle("The probability of the transformation of axions to X-rays in the magnet") +
    ggsave(&"out/transMagnetE_{year}.pdf")
  
  ############get the 1 and 2 sigma area ###################
  var pointR = pointdataR

  sort(pointR, system.cmp)

  var 
    sigma1 = (pointR.len.float * 0.68).round.int
    sigma2 = (pointR.len.float * 0.955).round.int
    sigmaAssign = newSeq[string](pointR.len)
  sigmaAssign.fill(0, sigma1 - 1, "sigma 1") 
  sigmaAssign.fill(sigma1, sigma2 - 1, "sigma 2")
  sigmaAssign.fill(sigma2, pointR.len - 1, "rest") 
  let rSigma1 = pointR[sigma1 - 1]
  let rSigma2 = pointR[sigma2 - 1]

  ############get the 1 and 2 sigma area before the Window###################

  var 
    sigma1Window = (axionsPassW.len.float * 0.68).round.int
    sigma2Window = (axionsPassW.len.float * 0.955).round.int
    sigmaAssignW = newSeq[string](pointR.len)
    rSigma1Window: float
    rSigma2Window: float
  if sigma2Window > pointR.len:
    sigmaAssignW.fill(0, sigma1Window - 1, "sigma 1") 
    sigmaAssignW.fill(sigma1Window, pointR.len - 1, "sigma 2")
    rSigma1Window = pointR[sigma1Window - 1]
  else:
    sigmaAssignW.fill(0, sigma1Window - 1, "sigma 1") 
    sigmaAssignW.fill(sigma1Window, sigma2Window - 1, "sigma 2")
    sigmaAssignW.fill(sigma2Window, pointR.len - 1, "rest") 
    rSigma1Window = pointR[sigma1Window - 1]
    rSigma2Window = pointR[sigma2Window - 1]


  #dfRad.mutate(fn {string -> string: "Sigma" ~ sigmaAssign[parseInt(`Idx`)]})
  let 
    dfRadOrg = dfRad.arrange("Radial component [mm]")
    pointdataRSig = dfRadOrg["Radial component [mm]"].toTensor(float)
    weightsSig = dfRadOrg["Transmission probability"].toTensor(float)
    pointDataXSig = dfRadOrg["x"].toTensor(float)
    pointDataYSig = dfRadOrg["y"].toTensor(float)
    sumWeights = sum(weightsSig)

  
  ##########same but for weighted values#####################
  var 
    sigma1Weight = (sumWeights * 0.68)
    sigma2Weight = (sumWeights * 0.955)
    sigmaAssignWeight = newSeq[string]((pointR.len.float * 0.63).round.int + 1)
    rSigma1W : float
    rSigma2W : float
    weightSum = sum(weightsSig[0..(pointR.len.float * 0.63).round.int])
  sigmaAssignWeight.fill(0, (pointR.len.float * 0.63).round.int, "sigma 1") 

  for i in (pointR.len.float * 0.63).round.int + 1..<pointR.len:
    weightSum += weightsSig[i]
    if weightSum < sigma1Weight:
      rSigma1W = pointR[i]
      sigmaAssignWeight.add("sigma 1") 
    elif weightSum < sigma2Weight and weightSum >= sigma1Weight:
      rSigma2W = pointR[i]
      sigmaAssignWeight.add("sigma 2")
    else:
      sigmaAssignWeight.add("rest") 
  echo rSigma1W, "vs ", rSigma1
  echo rSigma2W, "vs ", rSigma2


  let dfRadSig = seqsToDf({"Radial component [mm]": pointdataRSig,
                        "Transmission probability": weightsSig,
                        "x": pointDataXSig,
                       "y": pointDataYSig,
                       "Sigma" : sigmaAssignWeight,
                       "Sigma before window": sigmaAssignW})

  ggplot(dfRadSig, aes("Radial component [mm]", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma") +
    ggsave(&"out/radDistSig_{year}.pdf")

  ggplot(dfRadSig, aes("Radial component [mm]", fill = factor("Sigma before window"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma from before the window") +
    ggsave(&"out/radDistSigBeforeWindow_{year}.pdf")
  
  ggplot(dfRadSig, aes("x", "y", color = factor("Sigma"), weight = "Transmission probability")) +
    geompoint(size = some(0.5), alpha = some(0.1)) +
    ggtitle("X and Y") +
    ggsave(&"out/xy_{year}.pdf") 

  ggplot(dfRadSig, aes("y", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Y") +
    ggsave(&"out/y_{year}.pdf")


  #[let dfRadSigW = seqsToDf({"Radial component [mm]": pointdataRSig,
                        "Transmission probability": weightsSig,
                        "x": pointDataXSig,
                       "y": pointDataYSig,
                       "Sigma" : sigmaAssignWeight})

  ggplot(dfRadSigW, aes("Radial component [mm]", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma") +
    ggsave(&"out/radDistSigW_{year}.pdf")
  
  ggplot(dfRadSigW, aes("x", "y", fill = factor("Sigma"), weight = "Transmission probability")) +
    geompoint(size = some(0.5), alpha = some(0.1)) +
    ggtitle("X and Y") +
    ggsave(&"out/xyW_{year}.pdf") ]#

  #let fname2 = "extracted_from_aznar2015_llnl_telescope_eff_plot.csv"
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
    ggsave(&"out/FluxEafter_{year}.pdf")

  let dfFluxE3 = seqsToDf({"Axion energy [keV]": energiesPre})

  ggplot(dfFluxE3, aes("Axion energy [keV]")) +
    geom_histogram(binWidth = 0.1) +
    ylab("The flux before the experiment") +
    ggsave(&"out/FluxE_before_experiment_{year}.pdf")


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

  drawfancydiagrams("Axion Model Fluxfraction", heatmaptable2, 256, year, rSigma1W, rSigma2W) #rSigma1, rSigma2)

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
                        esCAST,
                         "2018",
                         "vacuum") # radiationCharacteristic = "axionRadiation::characteristic::sar"


