# stdlib
import std / [math, strutils, algorithm, random, sequtils, os, strformat, tables, sugar, strscans]

import ../axionMass/axionMassforMagnet

# nimble
import seqmath except linspace
import arraymancer except readCsv, linspace
import numericalnim, glm, ggplotnim, weave, cligen, unchained, parsetoml

##################rayTracer###############################

#degToRad(angle has to be done in Raytracer2014 for cos and sin

type
  ExperimentSetupKind = enum
    esCAST = "CAST"
    esBabyIAXO = "BabyIAXO"

  StageKind = enum
    skVacuum = "vacuum"
    skGas = "gas"

  WindowYearKind = enum
    wy2017 = "2017"
    wy2018 = "2018"
    wyIAXO = "BabyIAXO"

  CenterVectors = ref object
    centerEntranceCB: Vec3[float]
    centerExitCB: Vec3[float]
    centerExitPipeCBVT3: Vec3[float]
    centerExitPipeVT3XRT: Vec3[float]
    centerExitCBMagneticField: Vec3[float]
    centerEnterCBMagneticField: Vec3[float]
    centerXraySource: Vec3[float]
    centerCollimator: Vec3[float]
    centerSun: Vec3[float]

  ExperimentSetup* = ref object
    kind: ExperimentSetupKind
    stage: StageKind
    radiusCB*: mm
    RAYTRACER_LENGTH_COLDBORE*: mm
    RAYTRACER_LENGTH_COLDBORE_9T*: mm
    distXraySource*: mm
    radiusXraySource*: mm
    offAxXraySourceUp*: mm
    offAxXraySourceLeft*: mm
    lengthCol*: mm
    enXraySource*: keV
    activityXraySource*: GBq
    RAYTRACER_LENGTH_PIPE_CB_VT3*: mm
    radiusPipeCBVT3*: mm
    RAYTRACER_LENGTH_PIPE_VT3_XRT*: mm
    radiusPipeVT3XRT*: mm
    RAYTRACER_FOCAL_LENGTH_XRT*: mm
    distanceCBAxisXRTAxis*: mm
    RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW*: mm
    pipes_turned*: Degree
    optics_entrance: seq[MilliMeter]
    optics_exit: seq[MilliMeter]
    telescope_turned_x*: Degree
    telescope_turned_y*: Degree
    allThickness*: seq[MilliMeter]
    allR1*: seq[MilliMeter]
    allXsep*: seq[MilliMeter]
    allAngles*: seq[Degree]
    lMirror*: mm
    B*: T
    pGasRoom*: bar
    tGas*: K
    holeInOptics*: mm
    numberOfHoles: int
    holetype: string
    lateralDetector*: mm
    transversalDetector*: mm
    telescopeTransmission*: InterpolatorType[float] # should be `keV`, but cannot do that atm
    goldReflectivity*: Interpolator2DType[float]    # (angle / °, keV)

  DetectorSetupKind = enum
    dkInGrid2017 = "InGrid2017" # the setup as used in 2017
    dkInGrid2018 = "InGrid2018" # the setup as used in 2018
    dkInGridIAXO = "InGridIAXO"  # hypothetical setup for BabyIAXO (impoved window...)

  DetectorSetup = ref object # ref object for smaller size for parallel task buffer
    windowYear: WindowYearKind
    stripDistWindow: mm
    stripWidthWindow: mm
    detectorWindowAperture: mm # aperture of the detector window; isn't this just radius*2?
    theta: rad # rotation angle of the window (deduced from `windowYear`)
    radiusWindow: mm
    numberOfStrips: int
    openApertureRatio: float
    windowThickness: μm
    alThickness: μm
    depthDet: mm # depth (height) of the detector volume
    strongbackTransmission: InterpolatorType[float] # 1D interpolation for strongback transmission
    windowTransmission: InterpolatorType[float]     # 1D interpolation for window transmission
    gasAbsorption: InterpolatorType[float]          # 1D interpolation for gas absorption

  MaterialKind = enum
    mkSi3N4 = "Si3N4"
    mkSi = "Si"
    mkAr = "Ar"

  Axion = object
    passed: bool # indicates whether axion reached the detector
    passedTillWindow: bool
    hitNickel: bool
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
    energiesAx: keV
    energiesAxAll: keV
    energiesAxWindow: keV
    kinds: MaterialKind
    kindsWindow: MaterialKind
    transProbWindow: float
    transProbArgon: float
    transProbDetector: float
    transProbMagnet: float
    deviationDet: float
    shellNumber: int
    energiesPre: keV
    emratesPre: float
    reflect:float

  ConfigFlags = enum
    cfIgnoreDetWindow,   ## use to ignore the detector window absorbtion
    cfIgnoreGasAbs,      ## use to ignore the gas transmission
    cfIgnoreGoldReflect, ## use to ignore gold reflectivity (use 1.0)
    cfIgnoreConvProb,    ## use to ignore axion conversion probability
    cfXrayTest

################################
# VARIABLES from rayTracer.h
## WARNING: cannot be `const` at the moment, due to Nim compiler bug with distinct types
defUnit(GeV⁻¹)
let
  RAYTRACER_DISTANCE_SUN_EARTH = 1.5e14.mm # #ok
  radiusSun = 6.9e11.mm                    # #ok
  numberOfPointsSun = 1_000_000            #100000 for statistics   #37734 for CAST if BabyIaxo 10 mio  #26500960 corresponding to 100_000 axions at CAST, doesnt work
  # 1000000 axions that reach the coldbore then are reached after an operating time of 2.789 \times 10^{-5}\,\si{\second} for CAST


  roomTemp = 293.15.K
  mAxion = 0.0853#0.26978249412621896 #eV, corresponds to set p and T gas valus #0.4 #eV for example
  g_aγ = 2e-12.GeV⁻¹

## Chipregions#####

let
  CHIPREGIONS_CHIP_X_MIN = 0.0.mm
  CHIPREGIONS_CHIP_X_MAX = 66.0.mm #14.0.mm
  CHIPREGIONS_CHIP_Y_MIN = 0.0.mm
  CHIPREGIONS_CHIP_Y_MAX = 66.0.mm #14.0.mm
  CHIPREGIONS_CHIP_CENTER_X = CHIPREGIONS_CHIP_X_MAX / 2.0 #7.0.mm
  CHIPREGIONS_CHIP_CENTER_Y = CHIPREGIONS_CHIP_Y_MAX / 2.0 #7.0.mm
  CHIPREGIONS_GOLD_X_MIN = 4.5.mm
  CHIPREGIONS_GOLD_X_MAX = 9.5.mm
  CHIPREGIONS_GOLD_Y_MIN = 4.5.mm
  CHIPREGIONS_GOLD_Y_MAX = 9.5.mm
  CHIPREGIONS_SILVER_RADIUS_MAX = 4.5.mm
  CHIPREGIONS_BRONZE_RADIUS_MAX = 5.5.mm

################################

randomize(299792458)

proc initCenterVectors(expSetup: ExperimentSetup): CenterVectors =
  ## Initializes all the center vectors
  var centerSun = vec3(0.0)
  centerSun[0] = 0
  centerSun[1] = - (0.0 * 1.33e10)   ## first number number of millimeters at bore entrance
  centerSun[2] = - RAYTRACER_DISTANCE_SUN_EARTH.float

  var centerEntranceCB = vec3(0.0)
  centerEntranceCB[0] = 0
  centerEntranceCB[1] = -0.0
  centerEntranceCB[2] = 0 #coldboreBlockedLength # was 0 anyway

  var centerExitCBMagneticField = vec3(0.0)
  centerExitCBMagneticField[0] = 0
  centerExitCBMagneticField[1] = 0
  centerExitCBMagneticField[2] = expSetup.RAYTRACER_LENGTH_COLDBORE_9T.float

  var centerEnterCBMagneticField = vec3(0.0)

  var centerXraySource = vec3(0.0)
  centerXraySource[0] = expSetup.offAxXraySourceLeft.float
  centerXraySource[1] = expSetup.offAxXraySourceUp.float #250.0
  centerXraySource[2] = - (expSetup.distXraySource.float)

  var centerCollimator = vec3(0.0)
  centerCollimator[0] = expSetup.offAxXraySourceLeft.float
  centerCollimator[1] = expSetup.offAxXraySourceUp.float #250.0
  centerCollimator[2] = - (expSetup.distXraySource.float) + expSetup.lengthCol.float

  var centerExitCB = vec3(0.0)
  centerExitCB[0] = 0
  centerExitCB[1] = -0.0
  centerExitCB[2] = expSetup.RAYTRACER_LENGTH_COLDBORE.float

  var centerExitPipeCBVT3 = vec3(0.0)
  centerExitPipeCBVT3[0] = 0
  centerExitPipeCBVT3[1] = 0
  centerExitPipeCBVT3[2] = (expSetup.RAYTRACER_LENGTH_COLDBORE +
                            expSetup.RAYTRACER_LENGTH_PIPE_CB_VT3).float

  var centerExitPipeVT3XRT = vec3(0.0)
  centerExitPipeVT3XRT[0] = 0
  centerExitPipeVT3XRT[1] = 0
  centerExitPipeVT3XRT[2] = (expSetup.RAYTRACER_LENGTH_COLDBORE +
                             expSetup.RAYTRACER_LENGTH_PIPE_CB_VT3 +
                             expSetup.RAYTRACER_LENGTH_PIPE_VT3_XRT).float
  result = CenterVectors(centerEntranceCB: centerEntranceCB,
                         centerExitCB: centerExitCB,
                         centerExitPipeCBVT3: centerExitPipeCBVT3,
                         centerExitPipeVT3XRT: centerExitPipeVT3XRT,
                         centerExitCBMagneticField: centerExitCBMagneticField,
                         centerEnterCBMagneticField: centerEnterCBMagneticField,
                         centerXraySource: centerXraySource,
                         centerCollimator: centerCollimator,
                         centerSun: centerSun)

proc toRad(wyKind: WindowYearKind): float =
  ## Returns the radians corresponding to the angle of the detector window
  ## as it was installed in 2017 and 2018 of the CAST data taking campaign.
  ##
  ## Deduced from the calibration data & X-ray finger runs.
  ## TODO: Add reference to the sourcing of these numbers.
  case wyKind
  of wy2017:
    result = degToRad(10.8)
  of wy2018:
    result = degToRad(71.5)
  of wyIAXO:
    result = degToRad(20.0) # who knows

func conversionProb(B: Tesla, g_aγ: GeV⁻¹, length: MilliMeter): UnitLess =
  let L = length.mm.to(m)
  result = pow( g_aγ * B.toNaturalUnit() * L.toNaturalUnit() / 2.0, 2.0 )

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

proc `-`[T: SomeUnit](x: T): T =
  result = -1.0 * x

proc getRandomPointOnDisk(center: Vec3, radius: MilliMeter): Vec3 =
  ## This function gets a random point on a disk --> in this case this would
  ## be the exit of the coldbore ##
  var
    x = 0.0.mm
    y = 0.0.mm
    r = radius * sqrt(rand(1.0))
    angle = 360 * rand(1.0)
  x = cos(degToRad(angle)) * r
  y = sin(degToRad(angle)) * r
  result = vec3(x.float, y.float, 0.0) + center


proc getRandomPointFromSolarModel(center: Vec3, radius: MilliMeter,
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
  result = vec3(x.float, y.float, z.float) + center

template genGetRandomFromSolar(name, arg, typ, retTyp, body: untyped): untyped =
  ## This templates generates a procedure, which either returns a random energy
  ## for an event at a given radius or a random emission rate, biased
  ## by the emission rates at that radius.
  ## This only works if the energies to from are evenly distributed
  proc `name`(vectorInSun, center: Vec3, radius: MilliMeter,
              arg: typ,
              emRateCDFs: seq[seq[float]],
              ): retTyp =
    var
      rad = (vectorInSun - center).length.mm
      r = rad / radius # `UnitLess` radius
      iRad {.inject.}: int
      indexRad = (r - 0.0015) / 0.0005
    if indexRad - 0.5 > floor(indexRad):
      iRad = int(ceil(indexRad))
    else: iRad = int(floor(indexRad))
    # get the normalized (to 1) CDF for this radius
    let cdfEmRate = emRateCDFs[iRad]
    # sample an index based on this CDF
    let idx {.inject.} = cdfEmRate.lowerBound(rand(1.0))
    body

genGetRandomFromSolar(getRandomEnergyFromSolarModel,
                      energies, seq[keV],
                      keV):
  let energy = energies[idx]
  result = energy

genGetRandomFromSolar(getRandomEmissionRateFromSolarModel,
                      emissionRates, seq[seq[float]],
                      float):
  let emissionRate = emissionRates[iRad][idx]
  result = emissionRate

## The following are some functions to determine inetersection of the rays with the
## geometry (pipes, magnet, mirrors) of the setup ##
proc lineIntersectsCircle(point_1, point_2, center: Vec3,
                          radius: MilliMeter): bool =
  ## Now a function to see if the lines from the sun will actually intersect
  ## the circle area from the magnet entrance (called coldbore) etc. ##
  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  var intersect = vec3(0.0)
  intersect = (point_1 + lambda1 * vector) - center
  let r_xy_intersect = sqrt(intersect[0] * intersect[0] + intersect[1] *
      intersect[1])
  result = r_xy_intersect < radius.float

proc lineIntersectsObject(object_kind: string, point_1, point_2, center: Vec3,
                          radius: MilliMeter): bool =
  ## Now a function to see if the lines from the sun will actually intersect
  ## the circle area from the magnet entrance (called coldbore) etc. ##
  var
    vector = point_2 - point_1
    lambda1 = (center[2] - point_1[2]) / vector[2]
    intersect = (point_1 + lambda1 * vector) - center
    r_xy_intersect = sqrt(intersect[0] * intersect[0] + intersect[1] *
      intersect[1])
    intersect_turned = vec3(0.0)
  intersect_turned[0] = intersect[0] / sqrt(2.0) - intersect[1] / sqrt(2.0)
  intersect_turned[1] = intersect[0] / sqrt(2.0) + intersect[1] / sqrt(2.0)
  case object_kind
  of "circle":
    result = r_xy_intersect < radius.float
  of "cross":
    if (abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float * 16.0) or (abs(intersect[1]) < radius.float and abs(intersect[0]) < radius.float * 16.0):
      result = true
  of "star":
    if (abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float * 16.0) or (abs(intersect[1]) < radius.float and abs(intersect[0]) < radius.float * 16.0) or
    (abs(intersect_turned[0]) < radius.float and abs(intersect_turned[1]) < radius.float * 16.0) or (abs(intersect_turned[1]) < radius.float and abs(intersect_turned[0]) < radius.float * 16.0):
      result = true
  of "square":
    if abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float :
      result = true
  of "diamond":
    if abs(intersect_turned[0]) < radius.float and abs(intersect_turned[1]) < radius.float :
      result = true


proc getIntersectlineIntersectsCircle(point_1, point_2, center: Vec3): Vec3 =
  ## Now a function to get the intersection with one of the entrance cross sections
  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  result = point_1 + lambda1 * vector

proc lineIntersectsCylinderOnce(point_1: Vec3, point_2: Vec3, centerBegin: Vec3,
                                centerEnd: Vec3, radius: MilliMeter): bool =
  ## Also a function to know if the line intersected at least the whole magnet,
  ## and then only once, because else the axions would have just flown through ##
  ##
  ## TODO: can this *please* be merged with the proc below somehow?
  var 
    alpha_x = arcsin((centerEnd[0] - centerBegin[0]) / abs(centerBegin[2] - centerEnd[2]))
    alpha_y = arcsin((centerEnd[1] - centerBegin[1]) / abs(centerBegin[2] - centerEnd[2]))
    p_1 = point_1
    p_2 = point_2
    offset_x = 0.0
    offset_y = 0.0
  
  if abs(centerEnd[0]) <= abs(centerBegin[0]):
    offset_x = centerEnd[0]
  else: offset_x = centerBegin[0]
  if abs(centerEnd[1]) <= abs(centerBegin[1]):
    offset_y = centerEnd[1]
  else: offset_y = centerBegin[1]
  #echo alpha_y.radToDeg, " ", offset_y
  p_1[0] -= offset_x
  p_1[1] -= offset_y
  p_1[0] = p_1[0] * cos(alpha_x) + p_1[2] * sin(alpha_x)
  p_1[2] = p_1[2] * cos(alpha_x) - p_1[0] * sin(alpha_x)
  p_1[1] = p_1[1] * cos(alpha_y) - p_1[2] * sin(alpha_y)
  p_1[2] = p_1[2] * cos(alpha_y) + p_1[1] * sin(alpha_y)

  p_2[0] -= offset_x
  p_2[1] -= offset_y
  p_2[0] = p_2[0] * cos(alpha_x) + p_2[2] * sin(alpha_x)
  p_2[2] = p_2[2] * cos(alpha_x) - p_2[0] * sin(alpha_x)
  p_2[1] = p_2[1] * cos(alpha_y) - p_2[2] * sin(alpha_y)
  p_2[2] = p_2[2] * cos(alpha_y) + p_2[1] * sin(alpha_y)
  #echo alpha_y," ", point_2, " after ", p_2
  let
    vector = p_2 - p_1
    lambda_dummy = (-1000.0 - p_1[2]) / vector[2]
    dummy = p_1 + lambda_dummy * vector
    vector_dummy = p_2 - dummy
    factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
    p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
    q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - (radius*radius).float) / factor
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

proc getIntersectLineIntersectsCylinderOnce(
  point_1: Vec3, point_2: Vec3,
  centerBegin: Vec3, centerEnd: Vec3, radius: MilliMeter
     ): Vec3 =

  var 
    alpha_x = arcsin((centerEnd[0] - centerBegin[0]) / abs(centerBegin[2] - centerEnd[2]))
    alpha_y = arcsin((centerEnd[1] - centerBegin[1]) / abs(centerBegin[2] - centerEnd[2]))
    p_1 = point_1
    p_2 = point_2
    offset_x = 0.0
    offset_y = 0.0
  if abs(centerEnd[0]) <= abs(centerBegin[0]):
    offset_x = centerEnd[0]
  else: offset_x = centerBegin[0]
  if abs(centerEnd[1]) <= abs(centerBegin[1]):
    offset_y = centerEnd[1]
  else: offset_y = centerBegin[1]
  p_1[0] -= offset_x
  p_1[1] -= offset_y
  p_1[0] = p_1[0] * cos(alpha_x) + p_1[2] * sin(alpha_x)
  p_1[2] = p_1[2] * cos(alpha_x) - p_1[0] * sin(alpha_x)
  p_1[1] = p_1[1] * cos(alpha_y) - p_1[2] * sin(alpha_y)
  p_1[2] = p_1[2] * cos(alpha_y) + p_1[1] * sin(alpha_y)

  p_2[0] -= offset_x
  p_2[1] -= offset_y
  p_2[0] = p_2[0] * cos(alpha_x) + p_2[2] * sin(alpha_x)
  p_2[2] = p_2[2] * cos(alpha_x) - p_2[0] * sin(alpha_x)
  p_2[1] = p_2[1] * cos(alpha_y) - p_2[2] * sin(alpha_y)
  p_2[2] = p_2[2] * cos(alpha_y) + p_2[1] * sin(alpha_y)
  let
    vector = p_2 - p_1
    lambda_dummy = (-1000.0 - p_1[2]) / vector[2]
    dummy = p_1 + lambda_dummy * vector
    vector_dummy = p_2 - dummy
    factor = (vector_dummy[0]*vector_dummy[0] + vector_dummy[1]*vector_dummy[1])
    p = 2.0 * (dummy[0] * vector_dummy[0] + dummy[1]*vector_dummy[1]) / factor
    q = (dummy[0]*dummy[0] + dummy[1]*dummy[1] - (radius*radius).float) / factor
    lambda_1 = -p/2.0 + sqrt(p*p/4.0 - q)
    lambda_2 = -p/2.0 - sqrt(p*p/4.0 - q)
  var
    intersect_1 = dummy + lambda_1 * vector_dummy
    intersect_2 = dummy + lambda_2 * vector_dummy
    intersect_1_valid = (intersect_1[2] > centerBegin[2]) and
                        (intersect_1[2] < centerEnd[2])
    intersect_2_valid = (intersect_2[2] > centerBegin[2]) and
                        (intersect_2[2] < centerEnd[2])
  intersect_1[0] += offset_x
  intersect_1[1] += offset_y
  intersect_1[0] = intersect_1[0] * cos(-alpha_x) + intersect_1[2] * sin(-alpha_x)
  intersect_1[2] = intersect_1[2] * cos(-alpha_x) - intersect_1[0] * sin(-alpha_x)
  intersect_1[1] = intersect_1[1] * cos(-alpha_y) - intersect_1[2] * sin(-alpha_y)
  intersect_1[2] = intersect_1[2] * cos(-alpha_y) + intersect_1[1] * sin(-alpha_y)

  intersect_2[0] += offset_x
  intersect_2[1] += offset_y
  intersect_2[0] = intersect_2[0] * cos(-alpha_x) + intersect_2[2] * sin(-alpha_x)
  intersect_2[2] = intersect_2[2] * cos(-alpha_x) - intersect_2[0] * sin(-alpha_x)
  intersect_2[1] = intersect_2[1] * cos(-alpha_y) - intersect_2[2] * sin(-alpha_y)
  intersect_2[2] = intersect_2[2] * cos(-alpha_y) + intersect_2[1] * sin(-alpha_y)
  result = if (intersect_1_valid): intersect_1 else: intersect_2


proc getPixelValue(intersects: Vec3): Vec3 =
  const sizeViewfield = 48.0 #mm
  var intersectsPix = vec3(0.0)
  intersectsPix[0] = floor(intersects[0] / (sizeViewfield/1400.0)) + 700
  intersectsPix[1] = floor(intersects[1] / (sizeViewfield/1400.0)) + 700
  result = intersectsPix


## Some functions to include files from outside like the run file and the emissionrate/energy files ##

proc abs[T: SomeUnit](x: T): T =
  result = abs(x.float).T

proc findPosXRT*(pointXRT: Vec3, pointCB: Vec3,
                 r1, r2: MilliMeter, angle: Radian, lMirror, distMirr, uncer, sMin, sMax: MilliMeter): Vec3 =
  ## this is to find the position the ray hits the mirror shell of r1. it is after
  ## transforming the ray into a coordinate system, that has the middle of the
  ## beginning of the mirror cones as its origin
  var
    point = pointCB
    term: MilliMeter
    sMinHigh = sMin
    sMaxHigh = sMax
    pointMirror = vec3(0.0)
  let direc = pointXRT - pointCB
  var sValue = ((-direc[1] * point[0].mm * lMirror - direc[0] * point[1].mm * lMirror + direc[2].mm * r1.mm * sec(angle) - 
                direc[2].mm * r2 * sec(angle)).float - sqrt(pow((-direc[2].mm * r1.mm + direc[2].mm * r2 + 
                direc[1] * point[0].mm * lMirror * cos(angle) + direc[0] * point[1].mm * lMirror * cos(angle)).float, 2.0) - 
                4.0 * direc[0] * direc[1] * lMirror.float * cos(angle) * (distMirr * r1.mm - point[2].mm * r1.mm - 
                distMirr * r2 + point[2].mm * r2 + point[0].mm * point[1].mm * lMirror.float * cos(angle) + 
                lMirror * r1.mm * cos(angle)).float) * sec(angle))/(2.0 * direc[0] * direc[1] * lMirror.float)
  
  template calcVal(s: MilliMeter): untyped =
    ## Point + scalar * unit vector essentially. Hence no `mm` for direction.
    ## TODO: this should be handled differently...
    let res = sqrt((point[0].mm + s * direc[0]) * (point[0].mm + s * direc[0]) +
        (point[1].mm + s * direc[1]) * (point[1].mm + s * direc[1])) -
      ((r2 - r1) * (point[2].mm + s * direc[2] - distMirr) / (cos(angle) * lMirror))
    res
  var mid = (sMaxHigh + sMinHigh) / 2.0
  while abs(r1 - term) > 2.0 * uncer:
    if abs((sMinHigh - sMaxHigh)) < 1e-8.mm: break
    term = calcVal(mid)
    if abs(r1 - calcVal((sMinHigh + mid) / 2.0)) < abs(r1 - calcVal((sMaxHigh + mid) / 2.0)):
      # use lower half
      sMaxHigh = mid
      mid = (sMinHigh + mid) / 2.0
    else:
      # use upper half
      sMinHigh = mid
      mid = (sMaxHigh + mid) / 2.0
  pointMirror = point + mid.float * direc
  #echo sValue, " actual: ", mid.float
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
    distDet, xsepMiddle, dCBXray: MilliMeter,
    pipeAngle: Degree): Vec3 =

  ## To calculate the point in the detector window because the pipes are turned by 3 degree (angle here in rad)
  ## First switch into new coordinate system  with its origin in the middle of the telescope and z axis turned towards the detector
  let pipeRad = pipeAngle.to(Radian)
  var pointMirror2Turned = vec3(0.0)
  pointMirror2Turned[0] = (pointMirror2[0] * cos(pipeRad) +
                           pointMirror2[2] * sin(pipeRad)) - dCBXray.float
  pointMirror2Turned[1] = pointMirror2[1]
  pointMirror2Turned[2] = (pointMirror2[2] * cos(pipeRad) -
                           pointMirror2[0] * sin(pipeRad)) 
  var pointAfterMirror2Turned = vec3(0.0)
  pointAfterMirror2Turned[0] = (pointAfterMirror2[0] * cos(pipeRad) +
                                pointAfterMirror2[2] * sin(pipeRad)) - dCBXray.float
  pointAfterMirror2Turned[1] = pointAfterMirror2[1]
  pointAfterMirror2Turned[2] = (pointAfterMirror2[2] * cos(pipeRad) -
                                pointAfterMirror2[0] * sin(pipeRad)) 
  let vectorAfterMirror2 = pointAfterMirror2Turned - pointMirror2Turned
  ## Then the distance from the middle of the telescope to the detector can be calculated with the focal length
  ## Then n can be calculated as hown many times the vector has to be applied to arrive at the detector

  var distDet =  distDet / cos(pipeRad)
  var n = (distDet - pointMirror2Turned[2].mm) / vectorAfterMirror2[2].mm

  result = pointMirror2Turned + n.float * vectorAfterMirror2

## Now some functions for the graphs later, that store the data in heatmaps ##

proc prepareHeatmap(numberOfRows: int, numberOfColumns: int,
                    start_x: float, stop_x: float, start_y: float,
                    stop_y: float,
                    data_X: seq[float], data_Y: seq[float], weight1: seq[float],
                    norm: float64): Tensor[float] =
  ## This function prepares a heatmap out of given X and Y values with the z value
  ## (the number of entries in a certain pixel) as the weight of the event of the
  ## X and Y value
  # compute sizes based on number of coulmns / rows
  var stepsize_X = 0.0
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
  result = zeros[float]([numberOfRows, numberOfColumns])
  for i, value in data_X:
    var coord_X = floor((data_X[i] - start_x) / stepsize_X).int
    var coord_Y = floor((data_Y[i] - start_y) / stepsize_Y).int
    result[coord_Y, coord_X] = result[coord_Y, coord_X] + 1*weight1[i]/norm

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
  result = min(result, t.size - 1)

proc plotHeatmap(diagramtitle: string,
                 objectsToDraw: Tensor[float],
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
      zs[y * width + x] = objectsToDraw[y, x]

  #d.zmin = 0.0
  #d.zmax = 5e-22
  let offset = CHIPREGIONS_CHIP_CENTER_X.float
  let
    yr = linspace(- rSigma1, rSigma1, xs.len)
    yr2 = linspace(- rSigma2, rSigma2, xs.len)
    flux = zs.sum
  echo "The total flux arriving in the detector is: ", flux
  var df = seqsToDf({ "x" : xs,
                      "y" : ys,
                      "photon flux" : zs,
                      "yr0": yr,
                      "yr02": yr2})
    .mutate(f{float: "x-position [mm]" ~ `x` * CHIPREGIONS_CHIP_X_MAX.float / width.float},
            f{float: "y-position [mm]" ~ `y` * CHIPREGIONS_CHIP_Y_MAX.float / width.float},
            f{float: "xr" ~ sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + offset},
            f{float: "xrneg" ~ - sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + offset},
            f{float: "yr" ~ `yr0` + offset},
            f{float: "xr2" ~ sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + offset},
            f{float: "xrneg2" ~ - sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + offset},
            f{float: "yr2" ~ `yr02` + offset})
  template makeMinMax(knd, ax: untyped): untyped =
    template `knd ax`(): untyped =
      `CHIPREGIONS_GOLD ax knd` * width.float / CHIPREGIONS_CHIP_X_MAX.float
  makeMinMax(min, X)
  makeMinMax(max, X)
  makeMinMax(min, Y)
  makeMinMax(max, Y)

  let
    width = 720.0
    height = 586.0
  var customInferno = inferno()
  customInferno.name = "InfernoWithTransparent"
  customInferno.colors[0] = 0 shl 24 # transparent
  var customViridis = viridis()
  customViridis.name = "ViridisWithTransparent"
  customViridis.colors[0] = 0 shl 24 # transparent
  echo df
  #echo df.filter(f{`z` > 0.0})
  ggplot(df, aes("x-position [mm]", "y-position [mm]", fill = "photon flux")) +
    geom_raster() +
    scale_x_continuous() + scale_y_continuous() + scale_fill_continuous("photon flux") +
    scale_fill_gradient(customInferno) +
    #[geompoint(aes("xr", "yr"), color = some(parseHex("eab90c")), size = some(0.5)) +
    geompoint(aes("xrneg", "yr"), color = some(parseHex("eab90c")), size = some(0.5)) +
    geompoint(aes("xr2", "yr2"), color = some(parseHex("07529a")), size = some(0.5)) +
    geompoint(aes("xrneg2", "yr2"), color = some(parseHex("07529a")), size = some(0.5)) +]#
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = minY())) +
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = maxY())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = minX())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = maxX())) +
    #backgroundColor(parseHex("8cc7d4")) +
    #gridLineColor(parseHex("8cc7d4")) +
    #canvasColor(parseHex("8cc7d4")) +
    #theme_transparent() +
    margin(top = 2, left = 3, right = 6) +
    ggtitle(&"Simulated X-ray signal distribution on the detector chip with a total flux of {flux:.3e} events after 3 months") +
    ggsave(&"../out/axion_image_{year}.pdf", width = width, height = height)

  ggplot(df, aes("x-position [mm]", "y-position [mm]", fill = "photon flux")) +
    geom_raster() +
    scale_x_continuous() + scale_y_continuous() + scale_fill_continuous("photon flux") +
    #[geompoint(aes("xr", "yr"), color = some(parseHex("eab90c")), size = some(0.5)) +
    geompoint(aes("xrneg", "yr"), color = some(parseHex("eab90c")), size = some(0.5)) +
    geompoint(aes("xr2", "yr2"), color = some(parseHex("07529a")), size = some(0.5)) +
    geompoint(aes("xrneg2", "yr2"), color = some(parseHex("07529a")), size = some(0.5)) +]#
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = minY())) +
    # geom_line(aes = aes(xMin = minX(), xMax = maxX(), y = maxY())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = minX())) +
    # geom_line(aes = aes(yMin = minY(), yMax = maxY(), x = maxX())) +
    #annotate("sigma 1",                  # add our text annotation
          #x = 5.0, y = 0.1,                # at this location in 'data space'
          #backgroundColor = transparent) + # transparent background as we do manual TeX line breaks
    #canvasColor(parseHex("62bed3b0")) +
    theme_transparent() +
    ggtitle("Simulated X-ray signal distribution on the detector chip") +
    ggsave(&"../out/axion_image_{year}.png")

proc plotSolarModel(df: DataFrame) =
  ## A few plots for the solar model that are mainly for debugging.
  let df = df.mutate(f{"Energy" ~ (`Energy` * 9 + 1.0) * 0.001})
  ggplot(df, aes("Radius", "Flux")) +
    geom_point() +
    ggsave("../tmp/rad_flux.pdf")
  ggplot(df, aes("Energy", "Flux")) +
    geom_point() +
    ggsave("../tmp/energy_flux.pdf")
  ggplot(df, aes("Energy", "Flux", color = factor("Radius"))) +
    geom_line() +
    xlim(0, 0.1) +
    ggsave("../tmp/energy_rad_flux.pdf")
  let dfRad = df.group_by("Radius").summarize(f{float: "SumFlux" << sum(`Flux`)})
  echo dfRad
  ggplot(dfRad, aes("Radius", "SumFlux")) +
    geom_line() +
    ggsave("../tmp/flux_radius_sum.pdf")

############done with the functions, let's use them############
## NOTE: In principle it's a bit inefficient to re-parse the same config.toml file multiple times, but
## in the context of the whole ray tracing it doesn't matter. Makes the code a bit simpler.
const sourceDir = currentSourcePath().parentDir
proc parseResourcesPath(): string =
  ## parses the config.toml file containing the path to `resources` directory
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["resourcePath"].getStr

proc parseLlnlTelescopeFile(): string =
  ## parses the config.toml file containing the LLNL telescope efficiency filename
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["llnlEfficiency"].getStr

proc parseGoldFilePrefix(): string =
  ## parses the config.toml file containing the path to the gold reflectivity files
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["goldFilePrefix"].getStr

proc parseSolarModelFile(): string =
  ## parses the config.toml file containing the solar model filename
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["solarModelFile"].getStr

proc parseSetup(): (ExperimentSetupKind, DetectorSetupKind, StageKind) =
  ## parses the config.toml file containing the setup to compute the raytracing for
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result[0] = config["Setup"]["experimentSetup"].getStr.parseEnum[:ExperimentSetupKind]()
  result[1] = config["Setup"]["detectorSetup"].getStr.parseEnum[:DetectorSetupKind]()
  result[2] = config["Setup"]["stageSetup"].getStr.parseEnum[:StageKind]()

proc newExperimentSetup*(setup: ExperimentSetupKind,
                         stage: StageKind,
                         flags: set[ConfigFlags]): ExperimentSetup =
  # TODO: clean up, possibly make this into a toml file where one can
  # input different settings!
  case setup
  of esCAST:
    result = ExperimentSetup(
      kind: setup,
      stage: stage,
      radiusCB: 21.5.mm,
      RAYTRACER_LENGTH_COLDBORE: 9756.0.mm, # half B field to end of CB #ok
      RAYTRACER_LENGTH_COLDBORE_9T: 9260.0.mm, # half B field to half B field #ok
      distXraySource: 100.0.mm, #distance between the entrance of the magnet an a test Xray source
      radiusXraySource: 10.0.mm,
      offAxXraySourceUp: 200.0.mm,
      offAxXraySourceLeft: 0.0.mm,
      lengthCol: 50.0.mm,
      enXraySource: 1.0.keV,
      activityXraySource: 1.0.GBq,
      RAYTRACER_LENGTH_PIPE_CB_VT3: 2571.5.mm, # should stay the same #from beam pipe drawings #ok
      radiusPipeCBVT3: 39.64.mm, #30.0 # smallest aperture between end of CB and VT3
      RAYTRACER_LENGTH_PIPE_VT3_XRT: 150.0.mm, # from drawings #198.2 #mm from XRT drawing #ok
      radiusPipeVT3XRT: 35.0.mm, #25.0 # from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm #ok
      RAYTRACER_FOCAL_LENGTH_XRT: 1485.0.mm, #1300.0 # is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 # was the Telescope of 2014 (MPE XRT) also: Aperatur changed #ok
      distanceCBAxisXRTAxis: 0.0.mm, #62.1#58.44 # from XRT drawing #there is no difference in the axis even though the picture gets transfered 62,1mm down, but in the detector center
      RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW: 0.0.mm, # #no change, because don't know
      pipes_turned: 2.75.°, #degree # this is the angle by which the pipes before the detector were turned in comparison to the telescope
      optics_entrance: @[-83.0, 0.0, 0.0].mapIt(it.mm),
      optics_exit: @[-83.0, 0.0, 454.0].mapIt(it.mm),
      telescope_turned_x: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      telescope_turned_y: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
                             # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      allThickness: @[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
          0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2].mapIt(it.mm), ## the radii of the shells
      allR1: @[60.7095, 63.006, 65.606, 68.305, 71.105, 74.011, 77.027, 80.157,
          83.405, 86.775, 90.272, 93.902, 97.668, 101.576, 105.632].mapIt(it.mm), ## the radii of the shells
      allXsep: @[4.0, 4.171, 4.140, 4.221, 4.190, 4.228, 4.245, 4.288, 4.284,
          4.306, 4.324, 4.373, 4.387, 4.403, 4.481].mapIt(it.mm),
      #allR2: @[0.0, 60.731, 63.237, 65.838, 68.538, 71.339, 74.246, 77.263, 80.394, 83.642]
      allAngles: @[0.0, 0.579, 0.603, 0.628, 0.654, 0.680, 0.708, 0.737, 0.767,
          0.798, 0.830, 0.863, 0.898, 0.933, 0.970].mapIt(it.Degree), ## the angles of the mirror shells coresponding to the radii above
      lMirror: 225.0.mm, # Mirror length
      B: 9.0.T, # magnetic field of magnet
      pGasRoom: 1.0.bar, # pressure of the gas
      tGas: 1.7.K, #
      holeInOptics: 0.0.mm, #max 20.9.mm
      numberOfHoles: 5,
      holetype: "cross", #the type or shape of the hole in the middle of the optics
      lateralDetector: 0.0.mm, #lateral ofset of the detector in repect to the beamline
      transversalDetector: 0.0.mm #transversal ofset of the detector in repect to the beamline #0.0.mm #
    )
  of esBabyIAXO:
    result = ExperimentSetup(
      kind: setup,
      stage: stage,
      radiusCB: 500.0.mm, #350.0.mm,
                             # Change:
      RAYTRACER_LENGTH_COLDBORE: 11300.0.mm, # not sure if this is true but this is how its written on page 61 of the 2021 BabyIAXO paper
      RAYTRACER_LENGTH_COLDBORE_9T: 11000.0.mm, # I know it's not 9T here should be the actual length of pipe with a stable magnetic field; can't be same length
      distXraySource: 2000.0.mm, #88700.0.mm, #distance between the entrance of the magnet an a test Xray source
      radiusXraySource: 350.0.mm,
      offAxXraySourceUp: 0.0.mm,
      offAxXraySourceLeft: 0.0.mm,#4.5.mm, #
      lengthCol: 0.0.mm,
      enXraySource: 0.021.keV,
      activityXraySource: 0.125.GBq, #proposed source thing by Thomas #1.0.GBq,
      RAYTRACER_LENGTH_PIPE_CB_VT3: 225.0.mm, #300.0.mm, # not determined
      radiusPipeCBVT3: 370.0.mm, #mm smallest aperture between end of CB and VT4 # no Idea, I just made it wider than the coldbore
      RAYTRACER_LENGTH_PIPE_VT3_XRT: 250.0.mm, #300.0.mm, # not determined
      radiusPipeVT3XRT: 370.0.mm, # irrelevant, large enough to not loose anything # no idea
      RAYTRACER_FOCAL_LENGTH_XRT: 7500.0.mm, # # one possibility, the other is 5050 mm
      distanceCBAxisXRTAxis: 0.0.mm,
      RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW: 0.0.mm, # #no change, because don't know #good idea
      pipes_turned: 0.0.°, # this is the angle by which the pipes before the detector were turned in comparison to the telescope
      optics_entrance: @[0.0, -0.0, 0.0].mapIt(it.mm),
      optics_exit: @[0.0, -0.0, 600.0].mapIt(it.mm),
      telescope_turned_x: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      telescope_turned_y: -0.05.°, #the angle by which the telescope is turned in respect to the magnet
                             # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      #allR3: @[151.61, 153.88, 156.17, 158.48, 160.82, 163.18, 165.57, 167.98, 170.42, 172.88, 175.37, 177.88, 180.42, 183.14, 185.89, 188.67, 191.48,
          #194.32, 197.19, 200.09, 203.02, 206.03, 209.07, 212.14, 215.24, 218.37, 221.54, 224.74, 227.97, 231.24, 234.54, 237.87, 241.24, 244.85,
          #248.5, 252.19, 255.92, 259.68, 263.48, 267.32, 271.2, 275.12, 279.08, 283.09, 287.14, 291.38, 295.72, 300.11, 304.54, 309.02, 313.54,
          #318.11, 322.73, 327.4, 332.12, 336.88, 341.69, 346.55] #these are the real values but R3
      allThickness: @[0.468, 0.475, 0.482, 0.490, 0.497, 0.504, 0.511, 0.519, 0.526, 0.534, 0.542, 0.549, 0.557, 0.566, 0.574, 0.583, 0.591, 0.600, 0.609, 0.618,
          0.627, 0.636, 0.646, 0.655, 0.665, 0.675, 0.684, 0.694, 0.704, 0.714, 0.724, 0.735, 0.745, 0.756, 0.768, 0.779, 0.790, 0.802, 0.814, 0.826, 0.838, 0.850,
          0.862, 0.874, 0.887, 0.900, 0.913, 0.927, 0.941, 0.955, 0.968, 0.983, 0.997, 1.011, 1.026, 1.041, 1.055, 1.070].mapIt(it.mm),
      allR1: @[153.126, 155.419, 157.731, 160.065, 162.428, 164.812, 167.225, 169.66, 172.124, 174.608, 177.123, 179.658, 182.224, 184.971, 187.749, 190.556,
          193.394, 196.263, 199.161, 202.09, 205.05, 208.09, 211.16, 214.261, 217.392, 220.553, 223.755, 226.987, 230.249, 233.552, 236.885, 240.248, 243.652,
          247.298, 250.984, 254.711, 258.478, 262.276, 266.114, 269.992, 273.911, 277.87, 281.869, 285.92, 290.01, 294.292, 298.676, 303.109, 307.584, 312.108,
          316.674, 321.289, 325.955, 330.672, 335.439, 340.246, 345.104, 350.013].mapIt(it.mm), ## the radii of the shells closest to the magnet, now correct
      allXsep: @[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].mapIt(it.mm),
      allAngles: @[0.29, 0.294, 0.298, 0.303, 0.307, 0.312, 0.316, 0.321, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.355, 0.36, 0.366, 0.371, 0.377, 0.382, 0.388,
          0.393, 0.399, 0.405, 0.411, 0.417, 0.423, 0.429, 0.435, 0.441, 0.448, 0.454, 0.461, 0.467, 0.474, 0.481, 0.489, 0.496, 0.503, 0.51, 0.518, 0.525, 0.533,
          0.54, 0.548, 0.556, 0.564, 0.573, 0.581, 0.59, 0.598, 0.607, 0.616, 0.625, 0.634, 0.643, 0.652, 0.661].mapIt(it.Degree), ## the angles of the mirror shells coresponding to the radii above, now correct
      lMirror: 300.0.mm, # Mirror length
      B: 2.0.T, # magnetic field of magnet # Rather 2-3 T, not entirely homogeneous
      pGasRoom: 1.0.bar, #, pressure of the gas #for example P = 14.3345 mbar (corresponds to 1 bar at room temperature).
      tGas: 100.0.K, #293.15, # only Gas in BabyIAXO
      holeInOptics: 0.2.mm, #max 20.9.mm
      numberOfHoles: 1,
      holetype: "none", #the type or shape of the hole in the middle of the optics
      lateralDetector: 0.0.mm, #(sin(0.0.degToRad) * 7500.0).mm, #lateral ofset of the detector in repect to the beamline #0.0.mm #
      transversalDetector: (sin(0.0.degToRad) * 7500.0).mm #-0.0.mm # ##transversal ofset of the detector in repect to the beamline #0.0.mm #
    )

  ## TODO: this needs to be moved out of this procedure

  # TODO: this could be generalized to other telescopes until the multi layer stuff is implemented
  let resources = parseResourcesPath()
  let dfLlnl = readCsv(resources / parseLlnlTelescopeFile())
    .mutate(f{"Transmission" ~ idx("EffectiveArea[cm²]") /
      # total eff area of telescope = 1438.338mm² = 14.38338cm²
      (14.38338 - (4.3 * 0.2))}) # the last thing are the mirror seperators
  # TODO: take this out once happy
  when false:
    ggplot(dfTab["LLNL_transEff"], aes("Energy[keV]", "Transmission")) +
      geom_line() + ggsave("/tmp/transmission_llnl.pdf")

  var goldInterp: Interpolator2DType[float]
  if cfIgnoreGoldReflect notin flags:
    let prefix = resources / parseGoldFilePrefix()
    let goldFiles = collect(newSeq):
      for f in walkFiles(prefix & "*degGold0.25microns.csv"):
        let (success, angle) = scanTuple(f.dup(removePrefix(prefix)), "$fdeg")
        if not success: raise newException(IOError, "Could not parse input gold file: " & $f)
        (f, angle)
    var goldReflect = newTensor[float](0)
    for i, (f, angle) in goldFiles:
      let df = readCsv(f, sep = ' ', header = "#")
      if goldReflect.size == 0: # means first iteration, don't know # elements in gold files
        goldReflect = newTensor[float]([goldFiles.len, df.len])
      goldReflect[i, _] = df["Reflectivity", float].unsqueeze(0)
    let minAngle = goldFiles.mapIt(it[1]).min
    let maxAngle = goldFiles.mapIt(it[1]).max
    goldInterp = newBilinearSpline(goldReflect,
                                   (minAngle, maxAngle),
                                   (0.0, 15.0)) # 0 to 15 keV energy

  result.telescopeTransmission = newLinear1D(dfLlnl["Energy[keV]", float].toRawSeq,
                                             dfLlnl["Transmission", float].toRawSeq)
  result.goldReflectivity = goldInterp

defUnit(MilliMeter²)
proc sqrt(x: MilliMeter²): MilliMeter =
  ## We don't have sqrt as operator in unchained yet that does this automatically
  ## and errors if units are not perfect squares.
  sqrt(x.float).mm

proc calcWindowVals(radiusWindow: MilliMeter,
                    numberOfStrips: int,
                    openApertureRatio: float): tuple[width: MilliMeter,
                                                      dist: MilliMeter] =
  let
    totalArea = π * radiusWindow * radiusWindow
    areaOfStrips = totalArea * (1.0 - openApertureRatio)
    #width and distance between strips per strip;
    #the width on both sides is a whole width
    #(Don't know what to do about the additional string width)
    #not important at high strip number but at low like CAST
    dAndwPerStrip = radiusWindow * 2.0 / (numberOfStrips.float + 1.0)
  var
    lengthStrip: mm
    lengthAllStrips: mm

  for i in 0..(numberOfStrips/2).round.int - 1:
    lengthStrip = sqrt(radiusWindow * radiusWindow -
      (i.float * dAndwPerStrip + 0.5 * dAndwPerStrip) *
        (i.float * dAndwPerStrip + 0.5 * dAndwPerStrip)
    ) * 2.0
    lengthAllStrips = lengthAllStrips + lengthStrip
    echo lengthStrip

  lengthAllStrips = lengthAllStrips * 2.0

  let
    widthStrips = areaOfStrips / lengthAllStrips
    distStrips = dAndwPerStrip - widthStrips
  echo widthStrips
  echo distStrips
  result = (width: widthStrips, dist: distStrips)

proc newDetectorSetup*(setup: DetectorSetupKind): DetectorSetup =
  result = new DetectorSetup
  case setup
  of dkInGrid2017:
    result.windowYear = wy2017
    result.radiusWindow = 7.0.mm
    result.numberOfStrips = 4
    result.openApertureRatio = 0.838
    result.windowThickness = 0.3.μm
    result.alThickness = 0.02.μm # +- 0.007
    result.depthDet = 30.0.mm
  of dkInGrid2018:
    result.windowYear = wy2018
    result.radiusWindow = 7.0.mm
    result.numberOfStrips = 4
    result.openApertureRatio = 0.838
    result.windowThickness = 0.3.μm
    result.alThickness = 0.02.μm # +- 0.007
    result.depthDet = 30.0.mm
  of dkInGridIAXO:
    result.windowYear = wyIAXO
    result.depthDet = 30.0.mm
    result.radiusWindow = 7.0.mm #4.0.mm
    result.numberOfStrips = 4 #20 #maybe baby
    result.openApertureRatio = 0.838 #0.95 #
    result.windowThickness = 0.3.μm #0.1.μm #microns #options are: 0.3, 0.15 and 0.1
    result.alThickness = 0.02.μm #0.01.μm
  result.detectorWindowAperture = CHIPREGIONS_CHIP_X_MAX
  let (width, dist) = calcWindowVals(result.radiusWindow,
                                     result.numberOfStrips,
                                     result.openApertureRatio)
  result.stripWidthWindow = width
  result.stripDistWindow = dist
  ## TODO: clean this up! config.toml file!
  let resources = parseResourcesPath()
  let siNfile = resources / &"Si3N4Density=3.44Thickness={result.windowThickness.float:.1f}microns.tsv" #
  let siFile = resources / "SiDensity=2.33Thickness=200.microns.tsv"
  let detectorFile = resources / & "transmission-argon-30mm-1050mbar-295K.tsv"  #250
  let alFile = resources / &"AlDensity=2.7Thickness={result.alThickness.float:.2f}microns.tsv" #
  let dfSi  = readCsv(siFile, sep = ' ')
  let dfSiN = readCsv(siNfile, sep = ' ')
  var dfDet = readCsv(detectorFile, sep = ' ')
  let dfAl  = readCsv(alFile, sep = ' ')

  ## TODO: this needs to be moved out of this procedure
  var dfSB = dfSi
    .rename(f{"siFile" <- "Transmission"})
  dfSB["alFile"] = dfAl["Transmission", float]
  dfSB = dfSB.mutate(f{"Transmission" ~ `siFile` * `alFile`},
                     f{"Energy[keV]" ~ idx("PhotonEnergy(eV)") / 1000.0})
  var dfWd = dfSiN
    .rename(f{"siNFile" <- "Transmission"})
  dfWd["alFile"] = dfAl["Transmission", float]
  dfWd = dfWd.mutate(f{"Transmission" ~ `siNFile` * `alFile`},
                     f{"Energy[keV]" ~ idx("PhotonEnergy(eV)") / 1000.0})
  dfDet = dfDet
    .mutate(f{"Energy[keV]" ~ idx("PhotonEnergy(eV)") / 1000.0},
            f{"Absorption" ~ 1.0 - `Transmission`})
  result.strongbackTransmission = newLinear1D(dfSB["Energy[keV]", float].toRawSeq,
                                              dfSB["Transmission", float].toRawSeq)
  result.windowTransmission = newLinear1D(dfWd["Energy[keV]", float].toRawSeq,
                                          dfWd["Transmission", float].toRawSeq)
  result.gasAbsorption = newLinear1D(dfDet["Energy[keV]", float].toRawSeq,
                                     dfDet["Absorption", float].toRawSeq)
  result.theta = toRad(result.windowYear) # theta angle between window strips and horizontal x axis

template eval(interp: InterpolatorType[float], energy: keV): untyped =
  interp.eval(energy.float)

proc traceAxion(res: var Axion,
                centerVecs: CenterVectors,
                expSetup: ExperimentSetup,
                detectorSetup: DetectorSetup,
                emRates: seq[seq[float]],
                emRatesRadiusCumSum: seq[float],
                emRateCDFs: seq[seq[float]],
                energies: seq[keV],
                flags: set[ConfigFlags]
               ) =
  ## Get a random point in the sun, biased by the emission rate, which is higher
  ## at smalller radii, so this will give more points in the center of the sun ##
  var pointInSun = getRandomPointFromSolarModel(centerVecs.centerSun, radiusSun, emRatesRadiusCumSum)
  var weight = 1.0
  ## Get a random point at the end of the coldbore of the magnet to take all axions into account that make it to this point no matter where they enter the magnet ##
  var pointExitCBMagneticField = getRandomPointOnDisk(
      centerVecs.centerExitCBMagneticField, expSetup.radiusCB) #350.0.mm)
  let pointXraySource = getRandomPointOnDisk(
      centerVecs.centerXraySource, expSetup.radiusXraySource)
  let energyXraySource = expSetup.enXraySource

  ## Get a random energy for the axion biased by the emission rate ##
  var energyAx = getRandomEnergyFromSolarModel(
    pointInSun, centerVecs.centerSun, radiusSun, energies, emRateCDFs)
  if cfXrayTest in flags:
    pointInSun = pointXraySource
    energyAx = energyXraySource

    var centerSpot = vec3(0.0)
    if abs(expSetup.offAxXraySourceUp.float) > 50.0:
      centerSpot[0] = expSetup.offAxXraySourceLeft.float
      centerSpot[1] = expSetup.offAxXraySourceUp.float
    centerSpot[2] = expSetup.RAYTRACER_LENGTH_COLDBORE_9T.float
    #var radiusProjection = expSetup.radiusXraySource * 2.0 * (- centerVecs.centerXraySource[2] - expSetup.lengthCol.float).mm / expSetup.lengthCol + expSetup.radiusXraySource
    var radiusSpot = expSetup.holeInOptics.float
    if expSetup.holetype == "cross" or expSetup.holetype == "star":
      radiusSpot *= (expSetup.numberOfHoles.float + 16.0)
    elif expSetup.holetype == "none":
      radiusSpot = expSetup.radiusXraySource.float / 100.0 #* 3.0 #expSetup.radiusCB.float / 6.0 #
    else:
      radiusSpot *= (expSetup.numberOfHoles.float + 3.0)
    pointExitCBMagneticField[0] = pointInSun[0] + rand(0.1) - 0.05 #for parallel light
    pointExitCBMagneticField[1] = pointInSun[1] + rand(0.1) - 0.05 #for parallel light
    #pointExitCBMagneticField = getRandomPointOnDisk(centerSpot, (radiusSpot).mm) # for more statistics with hole through optics
    if not lineIntersectsCircle(pointInSun, pointExitCBMagneticField, centerVecs.centerCollimator, expSetup.radiusXraySource):

      
      return
    var xraysThroughHole = PI * radiusSpot * radiusSpot /
      (4.0 * PI * (- centerVecs.centerXraySource[2] + centerVecs.centerExitPipeVT3XRT[2]).mm *
      (- centerVecs.centerXraySource[2] + centerVecs.centerExitPipeVT3XRT[2]).mm) * expSetup.activityXraySource
    var testTime = 1_000_000 / (xraysThroughHole * 3600.0.s * 24.0)
    #echo "Days Testing ", testTime, " with a ", expSetup.activityXraySource, " source"
  #let emissionRateAx = getRandomEmissionRateFromSolarModel(
  #  pointInSun, centerVecs.centerSun, radiusSun, emRates, emRateCDFs
  #)
  ## Throw away all the axions, that don't make it through the piping system and therefore exit the system at some point ##



  let intersectsEntranceCB = lineIntersectsCircle(pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB, expSetup.radiusCB)
  var intersectsCB = false
  var pointInSunInSun = pointInSun - centerVecs.centerSun #
  #echo energyAx.float, " ", sqrt(pointInSunInSun[0].float * pointInSunInSun[0].float + pointInSunInSun[1].float * pointInSunInSun[1].float + pointInSunInSun[2].float * pointInSunInSun[2].float) / 6.9e11, " ", energyAx.float  * energyAx.float * (pointInSunInSun[0].float * pointInSunInSun[0].float + pointInSunInSun[1].float * pointInSunInSun[1].float + pointInSunInSun[2].float * pointInSunInSun[2].float) / 6.9e11 / 6.9e11
  res.emratesPre = 1.0 #* energyAx.float * energyAx.float * (pointInSunInSun[0].float * pointInSunInSun[0].float + pointInSunInSun[1].float * pointInSunInSun[1].float + pointInSunInSun[2].float * pointInSunInSun[2].float) / 6.9e11 / 6.9e11
  res.energiesPre = energyAx
  if (not intersectsEntranceCB):
    
    intersectsCB = lineIntersectsCylinderOnce(pointInSun,
        pointExitCBMagneticField, centerVecs.centerEntranceCB,
        centerVecs.centerExitCB, expSetup.radiusCB)
  if (not intersectsEntranceCB and not intersectsCB): return

  var intersect = vec3(0.0) #isnt't changed for axions that hit the entrance of the coldbore because the z value is 0 then anyways
  if (not intersectsEntranceCB): #generates problems with the weight because the weight is multiplied with the difference of the leght of the path of the particle and the legth of the coldbore
    intersect = getIntersectLineIntersectsCylinderOnce(
      pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB,
      centerVecs.centerExitCB, expSetup.radiusCB
    ) #pointInSun + ((centerVecs.centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
    #echo "doesn't intersect magnet entrance: ", intersect
  else:
    intersect = getIntersectlineIntersectsCircle(
      pointInSun,
      pointExitCBMagneticField, centerVecs.centerEntranceCB
    )
    #echo "does intersect magnet entrance: ", intersect

  ##get the length of the path of the axion in the magnetic field to get the probability of conversion later
  let pathCB = (pointExitCBMagneticField - intersect).length.mm
  var pointExitCB = vec3(0.0)
  #[if (not lineIntersectsCircle(pointInSun, pointExitCBMagneticField,
      centerVecs.centerExitCB, expSetup.radiusCB)): 
        echo "start"
        echo "exit magnet", pointExitCBMagneticField
        return]#

  pointExitCB = pointInSun + ((centerVecs.centerExitCB[2] - pointInSun[2]) / (
      pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)

  var pointExitPipeCBVT3 = vec3(0.0)

  #[if (not lineIntersectsCircle(pointExitCBMagneticField, pointExitCB,
      centerVecs.centerExitPipeCBVT3, expSetup.radiusPipeCBVT3)):
        echo "start"
        echo "exit magnet", pointExitCBMagneticField, "exit cb", pointExitCB
        return]#

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
  
  #echo centerVecs.centerCollimator, " ", pointExitPipeVT3XRT
  ###################from the CB (coldbore(pipe in Magnet)) to the XRT (XrayTelescope)#######################
  var vectorXRT = vectorBeforeXRT
  vectorXRT[0] = vectorXRT[0] * cos(expSetup.telescope_turned_x.to(Radian)) + vectorXRT[2] * sin(expSetup.telescope_turned_x.to(Radian))
  vectorXRT[2] = vectorXRT[2] * cos(expSetup.telescope_turned_x.to(Radian)) - vectorXRT[0] * sin(expSetup.telescope_turned_x.to(Radian))

  vectorXRT[1] = vectorXRT[1] * cos(expSetup.telescope_turned_y.to(Radian)) - vectorXRT[2] * sin(expSetup.telescope_turned_y.to(Radian))
  vectorXRT[2] = vectorXRT[2] * cos(expSetup.telescope_turned_y.to(Radian)) + vectorXRT[1] * sin(expSetup.telescope_turned_y.to(Radian))

  pointExitCB[0] -= expSetup.optics_entrance[0].float
  pointExitCB[1] -=  expSetup.optics_entrance[1].float
  pointExitCB[2] -=  centerVecs.centerExitPipeVT3XRT[2]
  pointExitCB[0] = pointExitCB[0] * cos(expSetup.telescope_turned_x.to(Radian)) + pointExitCB[2] * sin(expSetup.telescope_turned_x.to(Radian))
  pointExitCB[2] = pointExitCB[2] * cos(expSetup.telescope_turned_x.to(Radian)) - pointExitCB[0] * sin(expSetup.telescope_turned_x.to(Radian))

  pointExitCB[1] = pointExitCB[1] * cos(expSetup.telescope_turned_y.to(Radian)) - pointExitCB[2] * sin(expSetup.telescope_turned_y.to(Radian))
  pointExitCB[2] = pointExitCB[2] * cos(expSetup.telescope_turned_y.to(Radian)) + pointExitCB[1] * sin(expSetup.telescope_turned_y.to(Radian))

  var factor = (0.0 - pointExitCB[2]) / vectorXRT[2]
  
  #echo "before ", pointExitCB, " ", vectorBeforeXRT, " after ", pointExitCBXRT, " ", vectorXRT


  var pointEntranceXRT = pointExitCB + factor * vectorXRT
  
  vectorBeforeXRT = vectorXRT
  ## Coordinate transform from cartesian to polar at the XRT entrance
  var
    vectorEntranceXRTCircular = vec3(0.0)
  let
    radius1 = sqrt((pointEntranceXRT[0].mm) *
      (pointEntranceXRT[0].mm) + (pointEntranceXRT[1].mm) * (pointEntranceXRT[1].mm))
    phi_radius = arctan2(-pointEntranceXRT[1], (pointEntranceXRT[
        0])) #arccos((pointEntranceXRT[1]+d) / radius1)
    alpha = arctan(radius1 / expSetup.RAYTRACER_FOCAL_LENGTH_XRT)
    phi_flat = radtoDeg(arccos(pointEntranceXRT[0].mm / radius1))
  vectorEntranceXRTCircular[0] = radius1.float
  vectorEntranceXRTCircular[1] = phi_radius #in rad
  vectorEntranceXRTCircular[2] = alpha #in rad
  
  ### 3D spider structure 1st draft ###
  
  var factorSpider = (-85.0 - pointExitCB[2]) / vectorXRT[2]
  var pointEntrancSpider = pointExitCB + factorSpider * vectorXRT
  let
    radius1Spider = sqrt((pointEntrancSpider[0].mm) *
      (pointEntrancSpider[0].mm) + (pointEntrancSpider[1].mm) * (pointEntrancSpider[1].mm))
    phi_flatSpider = radtoDeg(arccos(pointEntrancSpider[0].mm / radius1Spider))
  #echo phi_flat, " ", phi_flatSpider, " ", pointEntrancSpider

  ## there is a 2mm wide graphite block between each glass mirror, to seperate them
  ## in the middle of the X-ray telescope. Return if hit
  ## BabyIAXO return if X-ray its the spider structure
  case expSetup.kind
  of esCAST:
    ## graphite blockers of the LLNL optics
    if pointEntranceXRT[1] <= 1.0 and pointEntranceXRT[1] >= -1.0: return
  of esBabyIAXO:
    ## here we have a spider structure for the XMM telescope:
    if vectorEntranceXRTCircular[0] <= 64.7: #and vectorEntranceXRTCircular[0] > expSetup.holeInOptics.float * expSetup.numberOfHoles.float: #doesnt really matter because there are no mirrors in the middle and these axions dont reach the window anyways
      for l in -(expSetup.numberOfHoles - ceil(expSetup.numberOfHoles.float / 2.0).int)..(expSetup.numberOfHoles - ceil(expSetup.numberOfHoles.float / 2.0).int):
        var centerHole = vec3(0.0)
        #var typeHole: string
        if l != 0:
          if abs(l) %% 2 == 0:
            centerHole[1] += 2.0 * l.float * expSetup.holeInOptics.float
          elif abs(l) %% 2 != 0:
            centerHole[0] += 2.0 * (l + (l / abs(l))).float * expSetup.holeInOptics.float
          #if abs(l) < 3: typeHole = expSetup.holetype
          #else: typeHole = "square"
        #else: typeHole = "square"
        if lineIntersectsObject(expSetup.holetype, pointExitCB, pointEntranceXRT, centerHole, expSetup.holeInOptics):
          weight = 1.0
          break
        else: weight = 0.0
    elif vectorEntranceXRTCircular[0] < 151.6 and vectorEntranceXRTCircular[0] > (151.6 - 20.9): #hole max 20.9 mm doesnt really matter because there are no mirrors in the middle and these axions dont reach the window anyways
      return
    elif vectorEntranceXRTCircular[0] > 64.7:
      for i in 0..16:
        if ((phi_flat >= (-1.25 + 22.5 * i.float) and phi_flat <= (1.25 + 22.5 * i.float))) or 
           ((phi_flatSpider >= (-1.25 + 22.5 * i.float) and phi_flatSpider <= (1.25 + 22.5 * i.float))): #spider strips (actually wider for innermost but doesnn't matter because it doesnt reach the window anyways)
          return
    #TODO: inner spider structure that doesnt matter
  if weight == 0.0: return

  ## Calculate the way of the axion through the telescope by manually reflecting the ray on the two mirror layers and then ending up before the detector ##

  var
    dist: seq[MilliMeter]
    r1 = 0.0.mm
    r2 = 0.0.mm
    r3 = 0.0.mm
    r4 = 0.0.mm
    r5 = 0.0.mm
    beta = 0.0.Radian ## in degree
    xSep = 0.0.mm
    h: int
    testXray = false
    lengthTelescope = (expSetup.lMirror + 0.5 * xSep) * cos(expSetup.allAngles[0].to(Radian)) + (expSetup.lMirror + 0.5 * xSep) * cos(3.0 * expSetup.allAngles[0].to(Radian))
    centerEndXRT = vec3(0.0)
  centerEndXRT[0] = 0.0 #lengthTelescope.float * tan(expSetup.telescope_turned.to(Radian)) thats not corrrect but the thing is already turned
  centerEndXRT[2] = lengthTelescope.float


  if vectorEntranceXRTCircular[0].mm > expSetup.allR1[expSetup.allR1.len - 1]: return
  for j in 0..<expSetup.allR1.len:
    # get rid of where the X-rays hit the glass frontal
    if vectorEntranceXRTCircular[0].mm > expSetup.allR1[j] and
        vectorEntranceXRTCircular[0].mm < (expSetup.allR1[j] + expSetup.allThickness[j]):
      return
    if expSetup.allR1[j] - vectorEntranceXRTCircular[0].mm > 0.0.mm:
      dist.add(expSetup.allR1[j] - vectorEntranceXRTCircular[0].mm)

  if min(dist) > 100.0.mm:
    testXray = true
    h = 50
    for l in -(expSetup.numberOfHoles - ceil(expSetup.numberOfHoles.float / 2.0).int)..(expSetup.numberOfHoles - ceil(expSetup.numberOfHoles.float / 2.0).int):
      var centerHole = centerEndXRT
      #var typeHole: string
      if l != 0:
        if abs(l) %% 2 == 0:
          centerHole[1] += 2.0 * l.float * expSetup.holeInOptics.float
        if abs(l) %% 2 != 0:
          centerHole[0] += 2.0 * (l + (l / abs(l))).float * expSetup.holeInOptics.float
        #if abs(l) < 3: typeHole = expSetup.holetype
        #else: typeHole = "square"
      #else: typeHole = "square"
      if lineIntersectsObject(expSetup.holetype, pointExitCB, pointEntranceXRT, centerHole, expSetup.holeInOptics):
        weight = 1.0
        break
      else: weight = 0.0
  if weight == 0.0: return
  if testXray == false:
    for k in 0..<expSetup.allR1.len:
      if min(dist) == expSetup.allR1[k] - vectorEntranceXRTCircular[0].mm:
        h = k


        r1 = expSetup.allR1[h]
        beta = expSetup.allAngles[h].to(Radian)
        xSep = expSetup.allXsep[h]
        r2 = r1 - expSetup.lMirror * sin(beta)
        r3 = r2 - 0.5 * xSep * tan(beta)
        r4 = r3 - 0.5 * xSep * tan(3.0 * beta)
        r5 = r4 - expSetup.lMirror * sin(3.0 * beta)

  if testXray == false and r1 == expSetup.allR1[0]: return

  var pointEntranceXRTZylKart = vec3(0.0)
  pointEntranceXRTZylKart[0] = pointEntranceXRT[0]
  pointEntranceXRTZylKart[1] = pointEntranceXRT[1]
  pointEntranceXRTZylKart[2] = pointEntranceXRT[2] #- centerVecs.centerExitPipeVT3XRT[2]

  var pointExitCBZylKart = vec3(0.0)
  pointExitCBZylKart[0] = pointExitCB[0]
  pointExitCBZylKart[1] = pointExitCB[1]
  pointExitCBZylKart[2] = pointExitCB[2] #- centerVecs.centerExitPipeVT3XRT[2]

  let
    beta3 = 3.0 * beta
    distanceMirrors = cos(beta) * (xSep + expSetup.lMirror)
    pointMirror1 = findPosXRT(pointEntranceXRTZylKart, pointExitCBZylKart, r1,
                              r2, beta, expSetup.lMirror, 0.0.mm, 0.001.mm, 1.0.mm, 1.1.mm)


  let
    vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRTZylKart,
        pointExitCBZylKart, pointMirror1, beta, "vectorAfter")
    pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    pointLowerMirror = findPosXRT(pointAfterMirror1, pointMirror1, expSetup.allR1[h-1] + expSetup.allThickness[h-1],
                              expSetup.allR1[h-1] - expSetup.lMirror * sin(expSetup.allAngles[h-1].to(Radian)) + expSetup.allThickness[h-1],
                              expSetup.allAngles[h-1].to(Radian),
                              expSetup.lMirror, 0.0.mm, 0.01.mm, 0.0.mm, 2.5.mm)
    pointMirror2 = findPosXRT(pointAfterMirror1, pointMirror1, r4, r5, beta3,
                              expSetup.lMirror, distanceMirrors, 0.01.mm, 0.0.mm, 2.5.mm)
  #echo r1, " ", expSetup.allR1[h-1], " ", r2, " ", expSetup.allR1[h-1] - expSetup.lMirror * sin(expSetup.allAngles[h-1].to(Radian))
  

  if pointMirror2[0] == 0.0 and pointMirror2[1] == 0.0 and pointMirror2[2] ==
      0.0: return ## with more uncertainty, 10% of the 0.1% we loose here can be recovered, but it gets more uncertain
  let
    vectorAfterMirrors = getVectoraAfterMirror(pointAfterMirror1, pointMirror1,
        pointMirror2, beta3, "vectorAfter")
    pointAfterMirror2 = pointMirror2 + 200.0 * vectorAfterMirrors
  let
    angle1 = getVectoraAfterMirror(pointEntranceXRTZylKart,
    pointExitCBZylKart, pointMirror1, beta, "angle")
    angle2 = getVectoraAfterMirror(pointAfterMirror1, pointMirror1,
        pointMirror2, beta3, "angle")
    alpha1 = angle1[1].round(2)
    alpha2 = angle2[1].round(2)
  #echo (angle1[1].degToRad) , " ", (r1.float - (expSetup.allR1[h-1] + expSetup.allThickness[h-1]).float) / (expSetup.lMirror.float - pointMirror1[2]) #radToDeg(arcsin(r1.float - (expSetup.allR1[h-1] + expSetup.allThickness[h-1]).float) / (expSetup.lMirror.float - pointMirror1[2]))

  # getting rid of the X-rays that hit the shell below 
  if testXray == false and tan(angle1[1].degToRad) > (r1.float - (expSetup.allR1[h-1] + expSetup.allThickness[h-1]).float) / (expSetup.lMirror.float - pointMirror1[2]):
    #if not (pointMirror1[2] +  0.0001 >  pointLowerMirror[2] and pointMirror1[2] -  0.0001 <  pointLowerMirror[2]):
    echo "hit Nickel"
    res.hitNickel = true
    return
    #echo pointLowerMirror, " ", pointMirror1, " ", pointMirror2, " ", angle1[1], " ", angle2[1], " ", h
  if (pointMirror1[2] +  0.001 >  pointMirror2[2] and pointMirror1[2] -  0.001 <  pointMirror2[2]): return
  ############################################# Mirrors end #################################################
  ## now get the points in the focal / detector plane

  ## because the detector is tuned in regards to the coldbore because it follows the direction of the telescope, set the origin to the detector window and
  ## turn the coordinate syste, to the detector for CAST

  var
    pointDetectorWindow = vec3(0.0)
    pointEndDetector = vec3(0.0)
    pointRadialComponent: float
    n: UnitLess
    n3: UnitLess
    distDet = distanceMirrors - 0.5 * expSetup.allXsep[8] * cos(beta) +
        expSetup.RAYTRACER_FOCAL_LENGTH_XRT -
        expSetup.RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW #distanc of the detector from the entrance of the optics
    d = - expSetup.optics_entrance[0]
  pointDetectorWindow = getPointDetectorWindow(
    pointMirror2,
    pointAfterMirror2,
    distDet, expSetup.allXsep[8], d, expSetup.pipes_turned
  )
  pointEndDetector = getPointDetectorWindow(
    pointMirror2, pointAfterMirror2,
    (distDet + detectorSetup.depthDet),
    expSetup.allXsep[8], d, expSetup.pipes_turned
  )
  

  ## TODO: what is this? 2^{x - y} ??
  res.deviationDet = sqrt(
    pow(2.0,
        (pointEndDetector[0] - pointDetectorWindow[0])) +
    pow(2.0,
        (pointEndDetector[1] - pointDetectorWindow[1])))

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
    distancePipe = (pointDetectorWindow[2] - pointExitCBZylKart[2]).mm.to(m) # needs to be in meter
  ## this is the transformation probability of an axion into a photon, if an axion
  ## flying straight through the magnet had one of 100%, angular dependency of the primakoff effect
  var
    transmissionMagnet: float

  case expSetup.stage
  of skVacuum:
    let probConversionMagnet = if cfIgnoreConvProb notin flags:
                                 conversionProb(expSetup.B, g_aγ, pathCB)
                               else:
                                 1.0
    transmissionMagnet = cos(ya) * probConversionMagnet
  of skGas:
    let
      pGas = expSetup.pGasRoom / roomTemp * expSetup.tGas
      ## TODO: use `unchained` in `axionMass` and avoid manual float conversions / use `to`
      effPhotonMass = effPhotonMass2(
        pGas.float, pathCB.to(m).float,
        expSetup.radiusCB.to(m).float,
        expSetup.tGas.float
      )
      probConversionMagnetGas =
        if cfIgnoreConvProb notin flags:
          axionConversionProb2(
            mAxion, energyAx, pGas.float,
            expSetup.tGas.float, pathCB.to(m).float,
            expSetup.radiusCB.to(m).float,
            g_aγ.float, expSetup.B.float
          ) # for setup including gas: functions are in axionmass/axionMassforMagnet
        else:
          1.0 # else perfect transmission
      absorbtionXrays = intensitySuppression2(
        energyAx, pathCB.to(m).float,
        distancePipe.float,
        pGas.float,
        expSetup.tGas.float,
        roomTemp.float) #room temperature in K
    #echo "axion mass in [eV] ", mAxion, " effective photon mass in [eV] ", effPhotonMass
    transmissionMagnet = cos(ya) * probConversionMagnetGas * absorbtionXrays # for setup with gas

  res.transmissionMagnets = transmissionMagnet
  res.yawAngles = ya

  

  var reflect = 0.0
  case expSetup.kind
  of esCAST:
    let transmissionTel = expSetup.telescopeTransmission.eval(energyAx)
    weight = (transmissionTel *
              transmissionTelescopePitch * transmissionTelescopeYaw *
              transmissionMagnet) #transmission probabilities times axion emission rate times the flux fraction
  of esBabyIAXO:
    if cfIgnoreGoldReflect notin flags:
      ## TODO: This needs to be replaced by a 2D interpolation!
      {.gcsafe.}:
        let
          reflectionProb1 = expSetup.goldReflectivity.eval(alpha1, energyAx.float)
          reflectionProb2 = expSetup.goldReflectivity.eval(alpha2, energyAx.float)
        res.reflect = reflectionProb1 * reflectionProb2
        weight = reflectionProb1 * reflectionProb2 * transmissionMagnet#also yaw and pitch dependend
    else:
      # without gold reflection just use perfect reflectivity
      weight = transmissionMagnet#also yaw and pitch dependend

  if testXray != false:
    if cfXrayTest notin flags:
      weight = transmissionMagnet
    n = (distDet - pointExitCB[2].mm) / (pointEntranceXRT - pointExitCB)[2].mm
    pointDetectorWindow = pointExitCB + n.float * (pointEntranceXRT - pointExitCB)
  pointDetectorWindow[0] -= expSetup.lateralDetector.float
  pointDetectorWindow[1] -= expSetup.transversalDetector.float
  if weight != 0:
    res.passedTillWindow = true
  #echo h, " ", pointEntranceXRT, " ", pointDetectorWindow
  ##Detector window:##
  if cfIgnoreDetWindow notin flags and sqrt(
      pointDetectorWindow[0].mm * pointDetectorWindow[0].mm +
      pointDetectorWindow[1].mm * pointDetectorWindow[1].mm
      ) > detectorSetup.radiusWindow:
      return

  else:
    if abs(pointDetectorWindow[0].mm) > CHIPREGIONS_CHIP_CENTER_X or abs(pointDetectorWindow[1].mm) > CHIPREGIONS_CHIP_CENTER_Y:
      return

  var pointDetectorWindowTurned = vec3(0.0)
  let theta = detectorSetup.theta
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
  # TODO: assignment here of the different kinds is obviously broken. Instead of having
  # one kinds field + the others we should have some additional field or something #made two assignments and now it works
  ## TODO: transmission of window material etc. can also be modeled using ray tracing.
  ## probability that transmission happens at all!
  let
    stripDistWindow = detectorSetup.stripDistWindow
    stripWidthWindow = detectorSetup.stripWidthWindow
  var transWindow = 0.0
  for i in 0..(detectorSetup.numberOfStrips/2).round.int - 1:
    if abs(y).mm > (1.0 * i.float + 0.5) * stripDistWindow + i.float * stripWidthWindow and
      abs(y).mm < (1.0 * i.float + 0.5) * stripDistWindow + (i.float + 1.0) * stripWidthWindow:
      transWindow = detectorSetup.strongbackTransmission.eval(energyAx)
      res.transProbWindow = transWindow
      res.transProbDetector = transWindow
      res.energiesAxAll = energyAx
      res.energiesAxWindow = energyAx
      res.kinds = mkSi
      res.kindsWindow = mkSi
      break
    else:
      transWindow = detectorSetup.windowTransmission.eval(energyAx)
      res.transprobWindow = transWindow
      res.transProbDetector = transWindow
      res.energiesAxAll = energyAx
      res.energiesAxWindow = energyAx
      res.kinds = mkSi3N4
      res.kindsWindow = mkSi3N4
  if cfIgnoreDetWindow notin flags:
    weight *= transWindow

  ## Get the total probability that the Xray will be absorbed by the detector and therefore detected:
  let absGasDet = detectorSetup.gasAbsorption.eval(energyAx)
  if cfIgnoreGasAbs notin flags:
    weight *= absGasDet
  res.transProbArgon = absGasDet
  res.transProbDetector = absGasDet
  res.energiesAxAll = energyAx
  res.kinds = mkAr
  res.energiesAx = energyAx
  res.shellNumber = h

  ###detector COS has (0/0) at the bottom left corner of the chip
  pointRadialComponent = sqrt(pointDetectorWindow[0]*pointDetectorWindow[0]+pointDetectorWindow[1]*pointDetectorWindow[1])
  res.pointdataR = pointRadialComponent
  pointDetectorWindow[0] = - pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X.float # for the view from the detector to the sun
  pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y.float

  if cfXrayTest notin flags:
    case expSetup.kind
    of esCAST:
      weight *= 3.585e3 * 3600.0 * 1.5 * 90.0 #for CAST at 10 mio, three months
    of esBabyIAXO:
      weight *= 9.5e6 * 3600.0 * 12.0 * 90.0 #9.5e6 * 3600.0 * 12.0 * 90.0 for three months at 1 mio axions at BabyIAXO]#
  ## assign final axion position & weight
  res.pointdataX = pointDetectorWindow[0]
  res.pointdataY = pointDetectorWindow[1]
  res.weights = weight
  res.weightsAll = weight


  ## TODO: the following are not used for anything!
  let
    gold = ( (pointDetectorWindow[0] >= CHIPREGIONS_GOLD_X_MIN.float) and (
        pointDetectorWindow[0] <= CHIPREGIONS_GOLD_X_MAX.float) and (
        pointDetectorWindow[1] >= CHIPREGIONS_GOLD_Y_MIN.float) and (
        pointDetectorWindow[1] <= CHIPREGIONS_GOLD_Y_MAX.float))
    r_xy = sqrt( ( (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X.float) * (
        pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X.float)) + ( (
        pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y.float) * (
        pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y.float))).mm
    silver = (r_xy <= CHIPREGIONS_SILVER_RADIUS_MAX.mm) and not gold
    bronze = not gold and not silver and (r_xy <= CHIPREGIONS_BRONZE_RADIUS_MAX)
    withinWindow = r_xy < detectorSetup.detectorWindowAperture / 2
    detector = ( (pointDetectorWindow[0] >= CHIPREGIONS_CHIP_X_MIN.float) and (
        pointDetectorWindow[0] <= CHIPREGIONS_CHIP_X_MAX.float) and (
        pointDetectorWindow[1] >= CHIPREGIONS_CHIP_Y_MIN.float) and (
        pointDetectorWindow[1] <= CHIPREGIONS_CHIP_Y_MAX.float))
  #if (energyAx < 0.3.keV) or (energyAx > 1.9.keV and energyAx < 2.1.keV):
    #echo energyAx, " ", weight #(transWindow), " tel: ", reflect, " det: ", 1.0 - transDet
  ## finally set the `passed` field to indicate this axion went all the way
  if weight != 0:
    res.passed = true
  #if energyAx > 0.5.keV:
    #res.passed = false



proc traceAxionWrapper(axBuf: ptr UncheckedArray[Axion],
                       bufLen: int,
                       centerVecs: CenterVectors,
                       expSetup: ExperimentSetup,
                       detectorSetup: DetectorSetup,
                       emRates: seq[seq[float]],
                       emRatesRadiusCumSum: seq[float],
                       emRateCDFs: seq[seq[float]],
                       energies: seq[keV],
                       flags: set[ConfigFlags]
                       ) =
  echo "Starting weave!"
  parallelFor iSun in 0 ..< bufLen:
    captures: { axBuf, centerVecs, expSetup, emRates, emRatesRadiusCumSum, emRateCDFs,
                detectorSetup,
                energies,
                flags }
    axBuf[iSun].traceAxion(centerVecs,
                           expSetup,
                           detectorSetup,
                           emRates, emRatesRadiusCumSum, emRateCDFs,
                           energies,
                           flags)

proc generateResultPlots(axions: seq[Axion],
                         windowYear: WindowYearKind) =
  ## Creates all plots we want based on the raytracing result
  let axionsPass = axions.filterIt(it.passed)
  echo "Passed axions ", axionsPass.len
  let axionsPassW = axions.filterIt(it.passedTillWindow)
  echo "Passed axions until the Window ", axionsPassW.len
  let hitNickelLayer = axions.filterIt(it.hitNickel)
  echo "Number of X-rays hitting nickel: ", hitNickelLayer.len

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
  extractAll(emratesPre)
  extractAll(transProbDetector)
  extractAll(transprobWindow)
  extractAll(kinds)
  extractAll(kindsWindow)
  extractPass(reflect)
  echo "Extracted all data!"
  ################################################################################
  ################################################################################
  ################################################################################

  #[let dfFluxTel = seqsToDf({"Axion energy [keV]": energiesAx.mapIt(it.float),
                            "Prob": reflect})

  let data = energiesAx.mapIt(it.float) # assume this is your data to be binned
  #let weights = reflect # assume these are our weights
  let binWidth = 0.0000666 # some width
  let nBins = 15000 #(data.max - data.min) / binWidth
  let (counts, bins) = histogram(data, bins = nBins) # unweighted
  let (countsW, _) = histogram(data, weights = reflect, bins = nBins)
  #echo data
  #echo counts
  let df = seqsToDf({"bins" : bins[0 ..< bins.high], "counts": counts, "countsW" : countsW})
    .mutate(f{"countsNorm" ~ `countsW` / `counts`})
  #echo df
  ggplot(df, aes("bins", "countsNorm")) +
    geom_histogram(stat = "identity") + # as we have already binned data!
    xlab("Axion energy [keV]") +
    ylab("Reflectivity") +
    margin(top = 2) +
    ggtitle("The reflectivity of the XMM optics depending on the incoming X-ray energy") +
    ggsave(&"../out/custom_histo_10.pdf") #to be used only with --ignoreDetWindow --ignoreGasAbs --ignoreConvProb


  echo dfFluxTel
  ggplot(dfFluxTel, aes("Axion energy [keV]", weight = "Prob")) +
    geom_histogram(binWidth = 0.0000666, lineWidth= some(1.2)) +
    xlim(0.0, 15.0) +
    ylim(0.0, 1000.0) +
    ylab("The flux before the experiment") +
    ggsave(&"../out/TelProb.pdf")]#

  #[let dfTransProb = seqsToDf({ "Axion energy [keV]": energiesAxAll.mapIt(it.float),
                               "Transmission Probability": transProbDetector,
                               "type": kinds.mapIt($it),
                               "Axion energy window[keV]":energiesAxWindow.mapIt(it.float),
                               "Transmission Probability window": transprobWindow,
                               "type window":kindsWindow.mapIt($it),
                               "Flux after experiment": weightsAll })

  ggplot(dfTransProb.arrange("Axion energy [keV]")) +
    geom_line(aes("Axion energy [keV]", "Transmission Probability",
             color = "type")) +
    geom_line(aes("Axion energy [keV]", "Transmission Probability window",
             color = "type window")) +
    geom_histogram(aes("Axion energy [keV]", weight = "Flux after experiment"), binWidth = 0.01) +
    ggtitle("The transmission probability for different detector parts") +
    ggsave(&"../out/TransProb_{windowYear}.pdf")

  let dfTransProbAr = seqsToDf({ "Axion energy [keV]": energiesAx.mapIt(it.float),
                                 "Transmission Probability": transProbArgon })
  ggplot(dfTransProbAr.arrange("Axion energy [keV]"),
         aes("Axion energy [keV]", "Transmission Probability")) +
    geom_line() +
    ggtitle("The transmission probability for the detector gas") +
    ggsave(&"../out/TransProbAr_{windowYear}.pdf")

  let dfDet = seqsToDf({ "Deviation [mm]": deviationDet,
                         "Energies": energiesAx.mapIt(it.float),
                         "Shell": shellNumber })
    .filter(f{Value: isNull(df["Shell"][idx]).toBool == false})

  ggplot(dfDet, aes("Deviation [mm]")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Deviation of X-rays detector entrance to readout") +
    ggsave(&"../out/deviationDet_{windowYear}.pdf")

  ggplot(dfDet, aes("Deviation [mm]", fill = factor("Shell"))) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Deviation of X-rays - detector entrance to readout") +
    ggsave(&"../out/deviationDet_stacked_{windowYear}.pdf")

  ggplot(dfDet, aes("Energies", fill = factor("Shell"))) +
    ggridges("Shell", overlap = 1.8) +
    geom_histogram(binWidth = 0.1, position = "identity") +
    ggtitle("X-ray energy distributions at detector") +
    ggsave(&"../out/energies_by_shell_{windowYear}.pdf", height = 600)

  ggplot(dfDet, aes("Deviation [mm]", fill = factor("Shell"))) +
    ggridges("Shell", overlap = 1.8) +
    geom_histogram(binWidth = 0.001, position = "identity") +
    ggtitle("Deviation of X-rays - detector entrance to readout") +
    ggsave(&"../out/deviationDet_ridges_{windowYear}.pdf", height = 600)]#





  let dfRad = seqsToDf({"Radial component [mm]": pointdataR,
                        "Transmission probability": weights,
                        "x": pointDataX,
                       "y": pointDataY})

  #echo dfRad.arrange("Radial component [mm]")
  ggplot(dfRad, aes("Radial component [mm]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions") +
    ggsave(&"../out/radialDistribution_{windowYear}.pdf")

  let dfFluxE = seqsToDf({ "Axion energy [keV]": energiesAx.mapIt(it.float),
                           "Transmission probability": weights })

  #[ggplot(dfFluxE, aes("Axion energy [keV]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001, lineWidth= some(1.2)) +
    #backgroundColor(parseHex("8cc7d4")) +
    #gridLineColor(parseHex("8cc7d4")) +
    #canvasColor(parseHex("8cc7d4")) +
    #theme_transparent() +
    ylab("photon flux") +
    #ylim(0.0, 0.0001) +
    xlim(0.0, 15.0) +
    ggtitle("Simulated photon flux depending on the energy of the axion") +
    ggsave(&"../out/fluxAfter_{windowYear}.pdf")

  ggplot(dfFluxE, aes("Axion energy [keV]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001, lineWidth= some(1.2)) +
    #backgroundColor(parseHex("8cc7d4")) +
    #gridLineColor(parseHex("8cc7d4")) +
    #canvasColor(parseHex("8cc7d4")) +
    #theme_transparent() +
    ylab("photon flux") +
    #ylim(0.0, 0.0001) +
    xlim(0.0, 1.0) +
    ggtitle("Simulated photon flux depending on the energy of the axion") +
    ggsave(&"../out/fluxAfterZoomed_{windowYear}.pdf")

  ggplot(dfFluxE, aes("Axion energy [keV]", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.00001, lineWidth= some(1.2)) +
    theme_transparent() +
    ggtitle("Simulated photon flux depending on the energy of the axion") +
    ggsave(&"../out/fluxAfter_{windowYear}.png")]#
  #[

  let dfXY = seqsToDf({"x": pointDataX,
                       "y": pointDataY,
                       "Transmission probability": weights})
    .mutate(f{"R" ~ sqrt(`x` * `x` + `y` * `y`)})


  ggplot(dfXY, aes("x", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("X and Y") +
    ggsave(&"../out/x_{windowYear}.pdf")

  ggplot(dfXY, aes("R", weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("R") +
    ggsave(&"../out/R_{windowYear}.pdf")

  let dfMag = seqsToDf({"Transmission probability": transmissionMagnets,
                       "Angles between path and magnetic field": yawAngles,
                       "Axion energy[keV]": energiesAx}).arrange("Angles between path and magnetic field")


  ggplot(dfMag, aes("Angles between path and magnetic field", "Transmission probability")) +
    geom_point(size = some(0.5), alpha = some(0.1)) +
    ylim(3.1e-24, 3.16e-24) +
    ggtitle("The probability of the transformation of axions to X-rays in the magnet") +
    ggsave(&"../out/transMagnet_{windowYear}.pdf")

  ggplot(dfMag, aes("Axion energy[keV]", "Transmission probability")) +
    geom_point(size = some(0.5), alpha = some(0.1)) +
    #ylim(3.1e-24, 3.16e-24) +
    ggtitle("The probability of the transformation of axions to X-rays in the magnet") +
    ggsave(&"../out/transMagnetE_{windowYear}.pdf")]#

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
  if sigma1Window > pointR.len:
    sigmaAssignW.fill(0, pointR.len - 1, "sigma 2")
  elif sigma1Window < pointR.len and sigma2Window > pointR.len:
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
    pointdataRSig = dfRadOrg["Radial component [mm]", float]
    weightsSig = dfRadOrg["Transmission probability", float]
    pointDataXSig = dfRadOrg["x", float]
    pointDataYSig = dfRadOrg["y", float]
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

  #[ggplot(dfRadSig, aes("Radial component [mm]", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma") +
    ggsave(&"../out/radDistSig_{windowYear}.pdf")

  ggplot(dfRadSig, aes("Radial component [mm]", fill = factor("Sigma before window"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma from before the window") +
    ggsave(&"../out/radDistSigBeforeWindow_{windowYear}.pdf")

  ggplot(dfRadSig, aes("x", "y", color = factor("Sigma"), weight = "Transmission probability")) +
    geompoint(size = some(0.5), alpha = some(0.1)) +
    ggtitle("X and Y") +
    ggsave(&"../out/xy_{windowYear}.pdf")]#
  let dfRadFilter = dfRadSig.filter(f{`x` < (CHIPREGIONS_CHIP_CENTER_X.float + 0.05)}).filter(f{`x` > (CHIPREGIONS_CHIP_CENTER_X.float - 0.05)})
  echo dfRadFilter
  ggplot(dfRadFilter, aes("y", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    xlab("y-position [mm]") +
    ylab("flux") +
    ggtitle("Simulated X-ray signal distribution along the y-axis") +
    ggsave(&"../out/y_{windowYear}.pdf")
  

  #[let dfRadSigW = seqsToDf({"Radial component [mm]": pointdataRSig,
                        "Transmission probability": weightsSig,
                        "x": pointDataXSig,
                       "y": pointDataYSig,
                       "Sigma" : sigmaAssignWeight})

  ggplot(dfRadSigW, aes("Radial component [mm]", fill = factor("Sigma"), weight = "Transmission probability")) +
    geom_histogram(binWidth = 0.001) +
    ggtitle("Radial distribution of the axions with sigma") +
    ggsave(&"../out/radDistSigW_{windowYear}.pdf")

  ggplot(dfRadSigW, aes("x", "y", fill = factor("Sigma"), weight = "Transmission probability")) +
    geompoint(size = some(0.5), alpha = some(0.1)) +
    ggtitle("X and Y") +
    ggsave(&"../out/xyW_{windowYear}.pdf") ]#

  #let fname2 = "extracted_from_aznar2015_llnl_telescope_eff_plot.csv"
  #let dfEnergyEff = toDf(readCsv(fname2, sep = ','))
  #  .mutate(fn {"Energies [keV]" ~ `xVals` * 8.0 + 1.0},
  #          fn {"Effective Area [cm^2]" ~ `yVals` *
  #              8.5}) #.mutate(f{int -> float: "Radii" ~ `Radii`.float * 0.0005 + 0.0015})

  #ggplot(dfEnergyEff, aes("Energies [keV]", "Effective Area [cm^2]")) +
  #  geom_line() +
  #  ggtitle("The telescope energy efficiency") +
  #  ggsave(&"../out/EnergyEff_{windowYear}.pdf")





  #[

  let dfFluxE2 = seqsToDf({ "Axion energy [keV]": energiesAx.mapIt(it.float),
                            "Flux after experiment": weights })
  dfFluxE2.write_csv(&"axion_gae_1e13_gagamma_{g_aγ.float}_flux_after_exp_N_{numberOfPointsSun}.csv")
  ggplot(dfFluxE2, aes("Axion energy [keV]", weight = "Flux after experiment")) +
    geom_histogram(binWidth = 0.001, lineWidth= some(1.2)) +
    ylab("The flux after the experiment") +
    ggsave(&"../out/FluxEafter_{windowYear}.pdf")

  let dfFluxE3 = seqsToDf({"Axion energy [keV]": energiesPre.mapIt(it.float),
                            "Flux before the experiment": emratesPre})

  ggplot(dfFluxE3, aes("Axion energy [keV]", weight = "Flux before the experiment")) +
    geom_histogram(binWidth = 0.0000666, lineWidth= some(1.2)) +
    xlim(0.0, 15.0) +
    ylab("The flux before the experiment") +
    ggsave(&"../out/FluxE_before_experiment_{windowYear}.pdf")



  ]#
  


  echo "all plots done, now to heatmap!"
  ## get the heatmaps out of the sequences of data X and data Y, first for the amount of data in one pixel ##
  ## compared to the overall amount and then the data in one pixel compared to the maximal amount of data in any pixel ##
  var
    beginX = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endX = CHIPREGIONS_CHIP_X_MAX.float  #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endY = CHIPREGIONS_CHIP_Y_MAX.float  #- distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareHeatmap(3000, 3000, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights,
      numberOfPointsSun.float) #colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareHeatmap(256, 256, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights, 1.0)
  var heatmaptable3 = prepareHeatmap(3000, 3000, beginX, endX, beginY, endY,
                                     pointdataX, pointdataY, weights, heatmaptable2.max) # if change number of rows: has to be in the maxVal as well
 # echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
 # echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]

  plotHeatmap("Axion Model Fluxfraction", heatmaptable2, 256, $windowYear, rSigma1W, rSigma2W) #rSigma1, rSigma2)

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

proc calculateFluxFractions(setup: ExperimentSetupKind,
                            detectorSetup: DetectorSetupKind,
                            stage: StageKind,
                            flags: set[ConfigFlags]) =
  let expSetup = newExperimentSetup(setup, stage, flags)
  let energies = linspace(0.001, 15.0, 15000).mapIt(it.keV)

  var
    integralNormalisation = 0.0
    integralTotal = 0.0
    integralDetector = 0.0
    integralBronze = 0.0
    integralSilver = 0.0
    integralGold = 0.0

  ## TODO: make the code use tensor for the emission rates!
  let resources = parseResourcesPath()
  var emRatesDf = readCsv(resources / parseSolarModelFile())
    .rename(f{"Radius" <- "dimension_1"}, f{"Energy" <- "dimension_2"}, f{"Flux" <- "value"})

  let emRatesTensor = emRatesDf["Flux", float]
    .reshape([emRatesDf.filter(fn {`Radius` == 0}).len, emRatesDf.filter(
        fn {`Energy` == 0}).len])

  let emRates = emRatesTensor
    .toRawSeq
    .reshape2D([emRatesTensor.shape[1], emRatesTensor.shape[0]])
  var
    emRatesRadiusCumSum: seq[float] = newSeq[float](emRates.len)
    emRateCDFs: seq[seq[float]] = newSeq[seq[float]](emRates.len)
    diffRadiusSum = 0.0

  template toCdf(x: untyped): untyped =
    let baseline = x[0]
    let integral = x[^1]
    x.mapIt( (it - baseline) / (integral - baseline) )

  for iRad in 0 ..< emRates.len:
    # emRates is seq of radii of energies
    var diffFlux = emRates[iRad]
    var diffSum = 0.0
    var radiusCumSum = newSeq[float](energies.len)
    for iEnergy in 0 ..< diffFlux.len:
      diffFlux[iEnergy] = diffFlux[iEnergy] * pow(energies[iEnergy].float, 2.0) * pow(iRad.float * 0.0005 + 0.0015, 2.0)
      diffSum += diffFlux[iEnergy]
      radiusCumSum[iEnergy] = diffSum
    diffRadiusSum += diffSum
    emRatesRadiusCumSum[iRad] = diffRadiusSum
    emRateCDFs[iRad] = radiusCumSum.toCdf()
  emRatesRadiusCumSum = emRatesRadiusCumSum.toCdf()

  when false:
    echo emRates[0].len, " ", emRates[0].len
    doAssert emRates[0].len == 15000
    var
      emRatesRadiusCumSum = emRates.mapIt(it.sum).cumSum() #zip(toSeq(0 ..< emRates.len), emRates).mapIt(it[1].sum * pow(it[0].float * 0.0005 + 0.0015, 2.0)).cumSum()
      emratesRadiusSum = emRates.mapIt(it.sum)
      emratesEnergySum: seq[float]

    # normalize to one
    emRatesRadiusCumSum.applyIt(it / emRatesRadiusCumSum[^1])
    for e in 1..15000:
      var emratesE = 0.0
      for r in 0..1967:
        emratesE += emRates[r][e-1] * (r.float * 0.0005 + 0.0015) * (r.float * 0.0005 + 0.0015)
        #echo e.float * 0.001, " ", (r.float * 0.0005 + 0.0015), " ", (r.float * 0.0005 + 0.0015) * (r.float * 0.0005 + 0.0015) * e.float * e.float * 0.000001
      emratesEnergySum.add(emratesE * e.float * e.float * 0.000001)
    #echo emratesEnergySum

    ggplot(seqsToDf({"Energies" : linspace(0.001, 15.0, 15000), "emrates" : emratesEnergySum}),
        aes("Energies", "emrates")) +
    geom_line() +
    xlim(0.0, 15.0) +
    ggsave("../out/energy_emrates.pdf")

    ## Compute all normalized CDFs of the emission rates for each radius
    var emRateCDFs = newSeq[seq[float]]()
    for iRad, f in emRates:
      var cdf = f.cumSum() #toSeq(0 ..< energies.len).mapIt(f[it] * (energies[it].to(eV).float * energies[it].to(eV).float) / (2 * Pi * Pi)).cumSum()
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
      es.add energyAx.float
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
  let detectorSetup = newDetectorSetup(detectorSetup)

  let centerVecs = expSetup.initCenterVectors()

  ## In the following we will go over a number of points in the sun, whose location and
  ## energy will be biased by the emission rate and whose track will be through the CAST
  ## experimental setup from 2018 at VT3
  var axions = newSeq[Axion](numberOfPointsSun)
  var axBuf = cast[ptr UncheckedArray[Axion]](axions[0].addr)
  echo "start"
  #echo emRatesRadiusSum.len
  echo emRates.len
  init(Weave)
  traceAxionWrapper(axBuf, numberOfPointsSun,
                    centerVecs,
                    expSetup,
                    detectorSetup,
                    emRates, emRatesRadiusCumSum,
                    emRateCDFs,
                    energies,
                    flags)
  exit(Weave)

  # walk the axions and determine `integralTotal` and `integral*`
  #if(gold and withinWindow): integralGold = integralGold + weight
  #if(silver and withinWindow): integralSilver = integralSilver + weight
  #if(bronze and withinWindow): integralBronze = integralBronze + weight
  #if(detector and withinWindow): integralDetector = integralDetector + weight

  generateResultPlots(axions, detectorSetup.windowYear)

proc main(ignoreDetWindow = false, ignoreGasAbs = false,
          ignoreConvProb = false, ignoreGoldReflect = false, xrayTest = false) =
  # check if the `config.toml` file exists, otherwise recreate from the default
  if not fileExists("config.toml"):
    let cdata = readFile("config_default.toml")
    writeFile("config.toml", cdata)

  var coldboreBlockedLength: float64
  coldboreBlockedLength = 0.0

  var flags: set[ConfigFlags]
  if ignoreDetWindow:   flags.incl cfIgnoreDetWindow
  if ignoreGasAbs:      flags.incl cfIgnoreGasAbs
  if ignoreConvProb:    flags.incl cfIgnoreConvProb
  if ignoreGoldReflect: flags.incl cfIgnoreGoldReflect
  if xrayTest: flags.incl cfXrayTest
  echo "Flags: ", flags

  let (esKind, dkKind, skKind) = parseSetup()

  calculateFluxFractions(esKind,
                         dkKind,
                         skKind,
                         flags) # radiationCharacteristic = "axionRadiation::characteristic::sar"


when isMainModule:
  dispatch main
