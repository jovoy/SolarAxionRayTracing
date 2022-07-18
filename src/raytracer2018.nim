# stdlib
import std / [math, strutils, algorithm, random, sequtils, os, strformat, tables, sugar, strscans, options]

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

  HoleType = enum
    htNone = "none"
    htCross = "cross"
    htStar = "star"
    htCircle = "circle"
    htSquare = "square"
    htDiamond = "diamond"

  TelescopeKind = enum
    tkLLNL, ## the LLNL telescope as used at CAST based on the nustar optics
    tkXMM,  ## the XMM Newton optics as they may be used in BabyIAXO
    tkCustomBabyIAXO, ## A hybrid optics specifically built for BabyIAXO mixing
                      ## a nustar like inner core with an XMM like outer core
    tkAbrixas, ## the Abrixas telescope as used at CAST
    tkOther

  StageKind = enum
    skVacuum = "vacuum"
    skGas = "gas"

  WindowYearKind = enum
    wy2017 = "2017"
    wy2018 = "2018"
    wyIAXO = "BabyIAXO"

  CenterVectors = ref object
    entranceCB: Vec3[float]
    exitCB: Vec3[float]
    exitPipeCBVT3: Vec3[float]
    exitPipeVT3XRT: Vec3[float]
    exitCBMagneticField: Vec3[float]
    enterCBMagneticField: Vec3[float]
    xraySource: Vec3[float]
    collimator: Vec3[float]
    sun: Vec3[float]

  ReflectivityKind = enum
    rkEffectiveArea, ## compute transmission probability from effective area (only used if
                     ## no reflectivity can be computed from first principles using telescope
                     ## layout and layer coatings)
    rkSingleCoating, ## used for telescopes with a single kind of coating on all layers (XMM, Abrixas)
    rkMultiCoating   ## used for telescopes with multiple different coatings (LLNL, custom BabyIAXO)

  Reflectivity = object
    case kind: ReflectivityKind
    of rkEffectiveArea:
      telescopeTransmission: InterpolatorType[float] # should be `keV`, but cannot do that atm
    of rkSingleCoating:
      reflectivity: Interpolator2DType[float]  # (angle / °, keV)
    of rkMultiCoating:
      ## Currently this assumes that the coatings are based on different layers and that the
      ## coatings change from inner most to outer most layer without repeating the same layer.
      ## Therefore, we have a mapping from layers to the correct interpolator via a lookup
      ## sequence (via linear scan), which stores the "boundary layers", i.e. the number of the
      ## first layer of the next type.
      layers: seq[int] ## e.g. @[14, 36, 60] layers 0-13 type 0, 14-35 type 1, ...,
                       ## given layer 18: layers.lowerBound(18) gives 1, implying to use reflectivity
                       ## interpolator at index 1
      reflectivities: seq[Interpolator2DType[float]]

  Magnet = object
    lengthColdbore: mm ## Length of the cold bore of the magnet
    B: T               ## Magnetic field strength of the magnet (assumed homogeneous)
    lengthB: mm        ## Length in which `B` applies
    radiusCB: mm       ## Radius of the coldbore of the magnet
    pGasRoom: bar      ## Pressure of gas inside the bore at room temperature
    tGas: K            ## Temperature of a possible gas inside the bore

  Telescope = object
    kind: TelescopeKind
    optics_entrance: seq[MilliMeter]
    optics_exit: seq[MilliMeter]
    telescope_turned_x: Degree
    telescope_turned_y: Degree
    allThickness: seq[MilliMeter]
    allR1: seq[MilliMeter]
    allXsep: seq[MilliMeter]
    allAngles: seq[Degree]
    lMirror: mm
    holeInOptics: mm
    numberOfHoles: int
    holeType: HoleType
    reflectivity: Reflectivity

  ## Information about an X-ray source installed for testing somewhere in front of the telescope
  TestXraySource = object
    active: bool  ## Whether the source is active (i.e. do we sample from the Sun or the source?)
    energy: keV   ## The energy of the X-ray source
    # distXraySource
    distance: mm  ## Distance of the X-ray source from the readout
    # radiusXraySource
    radius: mm    ## Radius of the X-ray source
    # offAxXraySourceUp
    offAxisUp: mm
    # offAxXraySourceLeft
    offAxisLeft: mm
    # activityXraySource
    activity: GBq ## The activity in `GBq` of the source
    lengthCol: mm ## Length of a collimator in front of the source

  ## An object that represents a single pipe
  Pipe = object
    length: mm
    radius: mm

  ## Stores information about the geometry of the pipes used between the
  ## magnet (coldbore), a gate valve (VT3) and the XRT
  Pipes = object
    # RAYTRACER_LENGTH_PIPE_CB_VT3
    # pipe connecting the coldbore itself to the VT3 gate valve
    coldBoreToVT3: Pipe
    # the pipe connecting the VT3 gate valve to the telescope (XRT)
    vt3ToXRT: Pipe
    # RAYTRACER_LENGTH_PIPE_VT3_XRT*: mm
    # radiusPipeVT3XRT*: mm
    distanceCBAxisXRTAxis: mm ## ??? set to 0.
    pipesTurned: Degree

  ## Describes where exactly the detector is placed. In "ideal" conditions the
  ## detector is placed exactly in the focal plane. In practice it describes
  ## the distance of the raytracing "camera" from the XRT. Therefore moving the
  ## detector installation changes the size / position of the image.
  DetectorInstallation = object
    # RAYTRACER_FOCAL_LENGTH_XRT*: mm
    distanceDetectorXRT: mm ## Distance from the XRT to the readout plane, i.e. where we record the raytracing image!
    # RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW
    distanceWindowFocalPlane: mm ## Distance from the detector window to the focal plane (can be used to offset
                                  ## `distanceDetectorXRT`)
    lateralShift: mm
    transversalShift: mm

  ExperimentSetup* = ref object
    kind: ExperimentSetupKind
    stage: StageKind
    magnet*: Magnet
    telescope*: Telescope
    testSource*: TestXraySource
    pipes*: Pipes
    detectorInstall*: DetectorInstallation

  DetectorSetupKind = enum
    dkInGrid2017 = "InGrid2017" # the setup as used in 2017
    dkInGrid2018 = "InGrid2018" # the setup as used in 2018
    dkInGridIAXO = "InGridIAXO"  # hypothetical setup for BabyIAXO (impoved window...)

  ## Describes the detector itself. What kind of gas it uses, the kind of window etc.
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

  ## Most of these field names should be renamed
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
    transmissionMagnet: float
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
    cfIgnoreReflection,  ## use to ignore reflectivity (use 1.0) of telescope layer coating
    cfIgnoreConvProb,    ## use to ignore axion conversion probability
    cfXrayTest,          ## if given reads the TestXraySource config values and uses it
    cfReadMagnetConfig,  ## if given uses configuration of magnet from config file
    cfReadDetInstallConfig ## if given uses configuration of detector install. from config file

  FullRaytraceSetup = object
    centerVecs: CenterVectors
    expSetup: ExperimentSetup
    detectorSetup: DetectorSetup
    energies: seq[keV]
    fluxRadiusCDF: seq[float]     ## CDF of the integrated flux per radius
                                  ## (indices correspond to radii)
    diffFluxCDFs: seq[seq[float]] ## CDF of each radius' differential flux per energy
                                  ## (indices correspond to radii & energies)
    flags: set[ConfigFlags]

################################
# VARIABLES from rayTracer.h
## WARNING: cannot be `const` at the moment, due to Nim compiler bug with distinct types
defUnit(GeV⁻¹)
let
  DistanceSunEarth = 1.5e14.mm # #ok
  RadiusSun = 6.9e11.mm                    # #ok
  NumberOfPointsSun = 1_000_000            #100000 for statistics   #37734 for CAST if BabyIaxo 10 mio  #26500960 corresponding to 100_000 axions at CAST, doesnt work
  # 1000000 axions that reach the coldbore then are reached after an operating time of 2.789 \times 10^{-5}\,\si{\second} for CAST

  RoomTemp = 293.15.K
  mAxion = 0.0853#0.26978249412621896 #eV, corresponds to set p and T gas valus #0.4 #eV for example
  g_aγ = 2e-12.GeV⁻¹

## Chipregions#####

let
  ChipXMin     = 0.0.mm
  ChipXMax     = 14.0.mm #14.0.mm
  ChipYMin     = 0.0.mm
  ChipYMax     = 14.0.mm #14.0.mm
  ChipCenterX  = ChipXMax / 2.0 #7.0.mm
  ChipCenterY  = ChipYMax / 2.0 #7.0.mm
  GoldXMin     = 4.5.mm
  GoldXMax     = 9.5.mm
  GoldYMin     = 4.5.mm
  GoldYMax     = 9.5.mm
  SilverRadius = 4.5.mm
  BronzeRadius = 5.5.mm

################################

randomize(299792458)

proc initCenterVectors(expSetup: ExperimentSetup): CenterVectors =
  ## Initializes all the center vectors, that is the center position of specific objects
  ## part of the raytracing in the global coordinate system.
  # Position of center of the Sun
  let sun = vec3(0.0,
                 - (0.0 * 1.33e10),   ## first number number of millimeters at bore entrance
                 - DistanceSunEarth.float)
  # position of the entrance of the magnet cold bore
  let entranceCB = vec3(0.0, -0.0, 0.0) #coldboreBlockedLength # was 0 anyway
  # position of beginning of magnetic field
  ## XXX: why is this the same as the `entranceCB`? According to the types the CB is a bit longer than
  ## the actual part of the magnetic field so shouldn't they differ?
  let enterCBMagneticField = vec3(0.0, 0.0, 0.0)
  # Exit of the CB magnetic field
  let exitCBMagneticField = vec3(0.0, 0.0, expSetup.magnet.lengthB.float)
  # position of the exit of the magnet cold bore
  let exitCB = vec3(0.0, -0.0, expSetup.magnet.lengthColdbore.float)
  # exit of the pipe connecting the coldbore to VT3
  let exitPipeCBVT3 = vec3(0.0, 0.0,
                           (expSetup.magnet.lengthColdbore +
                            expSetup.pipes.coldBoreToVT3.length).float)
  # exit of the pipe connecting VT3 to the telescope
  let exitPipeVT3XRT = vec3(0.0, 0.0,
                            (expSetup.magnet.lengthColdbore +
                             expSetup.pipes.coldBoreToVT3.length +
                             expSetup.pipes.vt3ToXRT.length).float)
  # position of an optional X-ray source for testing
  let xraySource = vec3(expSetup.testSource.offAxisLeft.float,
                        expSetup.testSource.offAxisUp.float, #250.0
                        - (expSetup.testSource.distance.float))
  # position of the collimator of the X-ray test source
  let collimator = vec3(expSetup.testSource.offAxisLeft.float,
                        expSetup.testSource.offAxisUp.float, #250.0
                        - (expSetup.testSource.distance.float) + expSetup.testSource.lengthCol.float)
  result = CenterVectors(entranceCB: entranceCB,
                         exitCB: exitCB,
                         exitPipeCBVT3: exitPipeCBVT3,
                         exitPipeVT3XRT: exitPipeVT3XRT,
                         exitCBMagneticField: exitCBMagneticField,
                         enterCBMagneticField: enterCBMagneticField,
                         xraySource: xraySource,
                         collimator: collimator,
                         sun: sun)

proc toRad(wyKind: WindowYearKind): float =
  ## Returns the radians corresponding to the angle of the detector window
  ## as it was installed in 2017 and 2018 of the CAST data taking campaign.
  ##
  ## Deduced from the calibration data & X-ray finger runs.
  ## TODO: Add reference to the sourcing of these numbers.
  case wyKind
  of wy2017, wy2018:
    result = degToRad(30.0)
  of wyIAXO:
    result = degToRad(20.0) # who knows

proc translateZ(vec: Vec3[float], distance: MilliMeter): Vec3[float] =
  result = vec
  result[2] += distance.float

proc rotateInX(vector: Vec3, angle: Radian, rotateOffset: MilliMeter = 0.0.mm): Vec3 =
  ## Rotation of a vector in x direction aka around the y axis counterclockwise when angle is positive and the y axis points towards observer
  ## Or rotation of the coordinate system the vector is in clockwise
  let vector = vector.translateZ(-rotateOffset)
  result = vec3(vector[0] * cos(angle) + vector[2] * sin(angle),
                vector[1],
                vector[2] * cos(angle) - vector[0] * sin(angle))
  result = result.translateZ(rotateOffset)

proc rotateInY(vector: Vec3, angle: Radian, rotateOffset: MilliMeter = 0.0.mm): Vec3 =
  ## Rotation of a vector in y direction aka around the x axis counterclockwise when angle is positive and the x axis points towards observer
  ## Or rotation of the coordinate system the vector is in clockwise
  let vector = vector.translateZ(-rotateOffset)
  result = vec3(vector[0],
                vector[1] * cos(angle) - vector[2] * sin(angle),
                vector[2] * cos(angle) + vector[1] * sin(angle))
  result = result.translateZ(rotateOffset)

proc rotateAroundZ(vector: Vec3, angle: Radian): Vec3 =
  ## Rotation of a vector around the z axis counterclockwise when angle is positive and the z axis points away from the observer
  ## Or rotation of the coordinate system the vector is in clockwise
  result = vec3(vector[0] * cos(angle) + vector[1] * sin(angle),
                vector[1] * cos(angle) - vector[0] * sin(angle),
                vector[2])

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
                                  fluxRadiusCDF: seq[float]): Vec3 =
  ## This function gives the coordinates of a random point in the sun, biased
  ## by the emissionrates (which depend on the radius and the energy) ##
  ##
  ## `fluxRadiusCDF` is the normalized (to 1.0) cumulative sum of the total flux per
  ## radius of all radii of the solar model.
  let
    angle1 = 360 * rand(1.0)
    angle2 = 180 * rand(1.0)
    ## random number from 0 to 1 corresponding to possible solar radii.
    randEmRate = rand(1.0)
    rIdx = fluxRadiusCDF.lowerBound(randEmRate)
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
              diffFluxCDFs: seq[seq[float]],
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
    let cdfEmRate = diffFluxCDFs[iRad]
    # sample an index based on this CDF
    let idx {.inject.} = cdfEmRate.lowerBound(rand(1.0))
    body

genGetRandomFromSolar(getRandomEnergyFromSolarModel,
                      energies, seq[keV],
                      keV):
  let energy = energies[idx]
  result = max(0.03.keV, energy)

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

proc lineIntersectsObject(objectKind: HoleType, point_1, point_2, center: Vec3,
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
  case objectKind
  of htCircle:
    result = r_xy_intersect < radius.float
  of htCross:
    if (abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float * 16.0) or
       (abs(intersect[1]) < radius.float and abs(intersect[0]) < radius.float * 16.0):
      result = true
  of htStar:
    if (abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float * 16.0) or
       (abs(intersect[1]) < radius.float and abs(intersect[0]) < radius.float * 16.0) or
       (abs(intersect_turned[0]) < radius.float and abs(intersect_turned[1]) < radius.float * 16.0) or
       (abs(intersect_turned[1]) < radius.float and abs(intersect_turned[0]) < radius.float * 16.0):
      result = true
  of htSquare:
    if abs(intersect[0]) < radius.float and abs(intersect[1]) < radius.float :
      result = true
  of htDiamond:
    if abs(intersect_turned[0]) < radius.float and abs(intersect_turned[1]) < radius.float :
      result = true
  of htNone:
    result = false ## XXX: or should this be true? guess not, as it wasn't handled before, i.e. was false

proc getIntersectlineIntersectsCircle(point_1, point_2, center: Vec3): Vec3 =
  ## Now a function to get the intersection with one of the entrance cross sections
  var vector = vec3(0.0)
  vector = point_2 - point_1
  var lambda1 = (center[2] - point_1[2]) / vector[2]
  result = point_1 + lambda1 * vector

type
  VecResult = tuple[vector: Vec3[float], valid: bool]
proc lineIntersectsCylinder(
  point_1: Vec3, point_2: Vec3,
  centerBegin: Vec3, centerEnd: Vec3, radius: MilliMeter
     ): tuple[inter1: VecResult, inter2: VecResult] =
  ## Computes the intersections of the ray with the cylinder (magnet bore)
  ## and returns a tuple of two intersection

  var
    alpha_x = arcsin((centerEnd[0] - centerBegin[0]) / abs(centerBegin[2] - centerEnd[2]))
    alpha_y = arcsin((centerEnd[1] - centerBegin[1]) / abs(centerBegin[2] - centerEnd[2]))
    offset_x = 0.0
    offset_y = 0.0
  if abs(centerEnd[0]) <= abs(centerBegin[0]):
    offset_x = centerEnd[0]
  else: offset_x = centerBegin[0]
  if abs(centerEnd[1]) <= abs(centerBegin[1]):
    offset_y = centerEnd[1]
  else: offset_y = centerBegin[1]

  let
    p_1 = rotateInY(rotateInX(point_1, alpha_x), alpha_y) - vec3(offset_x, offset_y, 0.0)
    p_2 = rotateInY(rotateInX(point_2, alpha_x), alpha_y) - vec3(offset_x, offset_y, 0.0)

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

  # get intersections unrotated
  var
    intersect_1 = dummy + lambda_1 * vector_dummy
    intersect_2 = dummy + lambda_2 * vector_dummy
  # undo "magnet" rotation
  intersect_1 += vec3(offset_x, offset_y, 0.0)
  intersect_1 = rotateInY(rotateInX(intersect_1, -alpha_x), -alpha_y)
  intersect_2 += vec3(offset_x, offset_y, 0.0)
  intersect_2 = rotateInY(rotateInX(intersect_2, -alpha_x), -alpha_y)
  # check if intersection was hit
  let
    intersect_1_valid = (intersect_1[2] > centerBegin[2]) and
                        (intersect_1[2] < centerEnd[2])
    intersect_2_valid = (intersect_2[2] > centerBegin[2]) and
                        (intersect_2[2] < centerEnd[2])
  result = (inter1: (vector: intersect_1, valid: intersect_1_valid),
            inter2: (vector: intersect_2, valid: intersect_2_valid))


proc lineIntersectsCylinderOnce(point_1: Vec3, point_2: Vec3, centerBegin: Vec3,
                                centerEnd: Vec3, radius: MilliMeter): bool =
  ## Also a function to know if the line intersected at least the whole magnet,
  ## and then only once, because else the axions would have just flown through ##
  let (inter1, inter2) = lineIntersectsCylinder(
    point_1, point_2, centerBegin, centerEnd, radius
  )
  if (inter1.valid and inter2.valid) or
    (not inter1.valid and not inter2.valid):
    result = false
  elif inter1.valid:
    result = true
  else:
    result = true
  ## TODO: is there another case than `elif inter1.valid`? Or why is that not the same
  ## branch?

proc getIntersectLineIntersectsCylinderOnce(
  point_1: Vec3, point_2: Vec3,
  centerBegin: Vec3, centerEnd: Vec3, radius: MilliMeter
     ): Vec3[float] =
  ## Returns the correct vector related to cutting the cylinder (coldbore)
  let (inter1, inter2) = lineIntersectsCylinder(
    point_1, point_2, centerBegin, centerEnd, radius
  )
  result = if (inter1.valid): inter1.vector else: inter2.vector

proc getPixelValue(intersects: Vec3): Vec3 =
  const sizeViewfield = 48.0 #mm
  var intersectsPix = vec3(0.0)
  intersectsPix[0] = floor(intersects[0] / (sizeViewfield/1400.0)) + 700
  intersectsPix[1] = floor(intersects[1] / (sizeViewfield/1400.0)) + 700
  result = intersectsPix

proc abs[T: SomeUnit](x: T): T =
  result = abs(x.float).T

proc findPosXRT*(pointXRT: Vec3, pointCB: Vec3,
                 r1, r2: MilliMeter, angle: Radian, lMirror, distMirr, uncer, sMin, sMax: MilliMeter,
                 name = ""): Vec3 =
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
  # calculate the values to solve for s with the p-q-formular, where p=b/a and q=c/a
  let
    k = tan(angle) * tan(angle)
    a = direc[0] * direc[0] + direc[1] * direc[1] - k * direc[2] * direc[2]
    b = 2.0 * (point[0] * direc[0] + r1.float * tan(angle) * direc[2] - k * (point[2] - distMirr.float) * direc[2])
    halfb = b / 2.0
    c = point[0] * point[0] + point[1] * point[1] - r1.float * r1.float + 2.0 * r1.float * tan(angle) * (point[2] - distMirr.float) -
        k * (point[2] - distMirr.float) * (point[2] - distMirr.float)

  # find nearest root that lies in acceptable range
  var
    s : float
    root1 = (-half_b - sqrt(half_b * half_b - a * c)) / a
    root2 = (-half_b + sqrt(half_b * half_b - a * c)) / a
  if point[2] + root1 * direc[2] > distMirr.float and point[2] + root1 * direc[2] < distMirr.float + lMirror.float * cos(angle):
    s = root1
  elif point[2] + root2 * direc[2] > distMirr.float and point[2] + root2 * direc[2] < distMirr.float + lMirror.float * cos(angle):
    s = root2
  else:
    s = 0.0


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
  pointMirror = point + s * direc
  #echo mid.float, " actual: ", s ," point: ", point + root2 * direc, " root2: ", root2, " target ", name
  result = pointMirror

proc getVectoraAfterMirror*(pointXRT, pointCB, pointMirror: Vec3,
                            angle: float, pointOrAngle: string): Vec3 =
  ## this is to find the vector after the reflection on the respective mirror
  ##
  ## TODO: clarify wheether this returns angles in rad or in degree!
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
  else:
    ## TODO: instead of having string arguments to determine what to do, split into
    ## separate functions!
    ##
    ## While doing that, also change the code to *not* use a Vec3 for the return type
    ## of the `"angle"` case. We only use the first argument of the returned vector!
    doAssert false, "Invalid choice of return value!"

proc getPointDetectorWindow(pointMirror2: Vec3, pointAfterMirror2: Vec3,
    distDet, xsepMiddle, dCBXray: MilliMeter,
    pipeAngle: Degree): Vec3 =

  ## To calculate the point in the detector window because the pipes are turned by 3 degree (angle here in rad)
  ## First switch into new coordinate system  with its origin in the middle of the telescope and z axis turned towards the detector
  let
    pipeRad = pipeAngle.to(Radian)
    pointMirror2Turned = rotateInX(pointMirror2, pipeRad) - vec3(dCBXray.float, 0.0, 0.0)
    pointAfterMirror2Turned = rotateInX(pointAfterMirror2, pipeRad) - vec3(dCBXray.float, 0.0, 0.0)
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
                 rSigma2: float,
                 suffix = "") =
  ## this function draws a diagram out a given heatmap ##
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
  let offset = ChipCenterX.float
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
    .mutate(f{float: "x-position [mm]" ~ `x` * ChipXMax.float / width.float},
            f{float: "y-position [mm]" ~ `y` * ChipYMax.float / width.float},
            f{float: "xr" ~ sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + offset},
            f{float: "xrneg" ~ - sqrt(rSigma1 * rSigma1 - `yr0` * `yr0`) + offset},
            f{float: "yr" ~ `yr0` + offset},
            f{float: "xr2" ~ sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + offset},
            f{float: "xrneg2" ~ - sqrt(rSigma2 * rSigma2 - `yr02` * `yr02`) + offset},
            f{float: "yr2" ~ `yr02` + offset})
  template makeMinMax(knd, ax: untyped): untyped =
    template `knd ax`(): untyped =
      `Gold ax knd` * width.float / ChipXMax.float
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
    ggtitle(&"Simulated X-ray signal distribution on the detector chip with a total flux of {flux:.3e} events after 3 months {suffix}") +
    ggsave(&"../out/axion_image_{year}_{suffix}.pdf", width = width, height = height)

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
    ggsave(&"../out/axion_image_{year}_{suffix}.png")

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

proc parseLlnlReflectivityFile(): string =
  ## parses the config.toml file containing the LLNL reflectivity filename
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["llnlReflFile"].getStr

proc parseGoldReflectivityFile(): string =
  ## parses the config.toml file containing the gold reflectivity filename as a H5 file
  let config = parseToml.parseFile(sourceDir / "config.toml")
  result = config["Resources"]["goldReflFile"].getStr

proc parseGoldFilePrefix(): string =
  ## parses the config.toml file containing the path to the gold reflectivity files
  ## as text data files directly from henke.gov
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

proc maybeParseMagnetConfig(flags: set[ConfigFlags]): Option[Magnet] =
  ## parses the `Magnet` configuration from the `config.toml` file if
  ## the `cfReadMagnetConfig` flag is set or the `useConfig`
  ## from the `Magnet` config table is set to `true`.
  ##
  ## Returns `some[Magnet]` if parsing took place or `none[Magnet]` if
  ## we should use the experiment specific magnet.
  let cfg = parseToml.parseFile(sourceDir / "config.toml")["Magnet"]
  if cfReadMagnetConfig in flags or cfg["useConfig"].getBool:
    result = some(
      Magnet(B:              cfg["B"].getFloat.T,
             lengthB:        cfg["lengthB"].getFloat.mm,
             radiusCB:       cfg["radiusCB"].getFloat.mm,
             lengthColdbore: cfg["lengthColdbore"].getFloat.mm,
             pGasRoom:       cfg["pGasRoom"].getFloat.bar,
             tGas:           cfg["tGas"].getFloat.K
      )
    )
  else:
    result = none[Magnet]() # default, but let's be explicit

proc maybeParseTestXraySource(flags: set[ConfigFlags]): Option[TestXraySource] =
  ## parses the `TestXraySource` configuration from the `config.toml` file if
  ## the `cfReadTestXraySourceConfig` flag is set or the `useConfig`
  ## from the `TestXraySource` config table is set to `true`.
  ##
  ## Returns `some[TestXraySource]` if parsing took place or `none[TestXraySource]` if
  ## we should use the experiment specific magnet.
  let cfg = parseToml.parseFile(sourceDir / "config.toml")["TestXraySource"]
  if cfXrayTest in flags or cfg["useConfig"].getBool:
    result = some(
      TestXraySource(
        active:      cfg["active"].getBool,
        energy:      cfg["energy"].getFloat.keV,
        distance:    cfg["distance"].getFloat.mm,
        radius:      cfg["radius"].getFloat.mm,
        offAxisUp:   cfg["offAxisUp"].getFloat.mm,
        offAxisLeft: cfg["offAxisLeft"].getFloat.mm,
        activity:    cfg["activity"].getFloat.GBq,
        lengthCol:   cfg["lengthCol"].getFloat.mm
      )
    )
  else:
    result = none[TestXraySource]() # default, but let's be explicit

proc maybeParseDetectorInstallation(flags: set[ConfigFlags]): Option[DetectorInstallation] =
  ## parses the `DetectorInstallation` configuration from the `config.toml` file if
  ## the `cfReadDetInstallConfig` flag is set or the `useConfig`
  ## from the `DetectorInstallation` config table is set to `true`.
  ##
  ## Returns `some[DetectorInstallation]` if parsing took place or `none[DetectorInstallation]` if
  ## we should use the experiment specific magnet.
  let cfg = parseToml.parseFile(sourceDir / "config.toml")["DetectorInstallation"]
  if cfReadDetInstallConfig in flags or cfg["useConfig"].getBool:
    result = some(
      DetectorInstallation(
        distanceDetectorXRT: cfg["distanceDetectorXRT"].getFloat.mm,
        distanceWindowFocalPlane: cfg["distanceWindowFocalPlane"].getFloat.mm,
        lateralShift: cfg["lateralShift"].getFloat.mm,
        transversalShift: cfg["transversalShift"].getFloat.mm
      )
    )
  else:
    result = none[DetectorInstallation]() # default, but let's be explicit

proc initMagnet(setup: ExperimentSetupKind, flags: set[ConfigFlags]): Magnet =
  let magnetOpt = maybeParseMagnetConfig(flags)
  if magnetOpt.isSome:
    return magnetOpt.get
  case setup
  of esCAST:
    ## Build a CAST LHC dipole prototype magnet
    result = Magnet(
      B: 9.0.T, # magnetic field of magnet
      radiusCB: 21.5.mm, # radius of the cold bore
      lengthColdbore: 9756.0.mm, # half B field to end of CB #ok
      lengthB: 9260.0.mm, # half B field to half B field #ok
      pGasRoom: 1.0.bar, # pressure of the (optional) gas
      tGas: 1.7.K # temperature of (optional) gas
    )
  of esBabyIAXO:
    ## Build a BabyIAXO magnet
    result = Magnet(
      B: 2.0.T, # magnetic field of magnet # Rather 2-3 T, not entirely homogeneous
      radiusCB: 500.0.mm, #350.0.mm,
                             # Change:
      lengthColdbore: 11300.0.mm, # not sure if this is true but this is how its written on page 61 of the 2021 BabyIAXO paper
      lengthB: 11000.0.mm, # I know it's not 9T here should be the actual length of pipe with a stable magnetic field; can't be same length
      pGasRoom: 1.0.bar, #, pressure of the gas #for example P = 14.3345 mbar (corresponds to 1 bar at room temperature).
      tGas: 100.0.K #293.15, # only Gas in BabyIAXO
    )

proc initPipes(setup: ExperimentSetupKind): Pipes =
  case setup
  of esCAST:
    result = Pipes(
      ## the next variable doesn't make sense. It's not 2.57m from the cold bore to VT3!
      coldBoreToVT3: Pipe(length: 2571.5.mm, # should stay the same #from beam pipe drawings #ok
                          radius: 39.64.mm), #30.0 # smallest aperture between end of CB and VT3
      vt3ToXRT: Pipe(length: 150.0.mm, # from drawings #198.2 #mm from XRT drawing #ok
                     radius: 35.0.mm), #25.0 # from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm #ok
      pipesTurned: 2.75.°, #degree # this is the angle by which the pipes before the detector were turned in comparison to the telescope
      distanceCBAxisXRTAxis: 0.0.mm #62.1#58.44 # from XRT drawing #there is no difference in the axis even though the picture gets transfered 62,1mm down, but in the detector center
    )
  of esBabyIAXO:
    result = Pipes(
      coldBoreToVT3: Pipe(length: 225.0.mm, #300.0.mm, # not determined
                          radius: 370.0.mm), #mm smallest aperture between end of CB and VT4 # no Idea, I just made it wider than the coldbore
      vt3ToXRT: Pipe(length: 250.0.mm, #300.0.mm, # not determined
                     radius: 370.0.mm), # irrelevant, large enough to not loose anything # no idea
      pipesTurned: 0.0.°, # this is the angle by which the pipes before the detector were turned in comparison to the telescope
      distanceCBAxisXRTAxis: 0.0.mm
    )

import nimhdf5 except linspace # only needed for here
proc initReflectivity(optics: TelescopeKind): Reflectivity =
  ## Inits the `reflectivity` field of an `ExperimentSetup`, i.e. the reflectivities of
  ## the telescope layers / the effecitve reflectivity based on the effective area, if applicable
  case optics
  of tkLLNL:
    result = Reflectivity(
      kind: rkMultiCoating,
      layers: @[2, 2+3, 2+3+4, 2+3+4+5] # layers of LLNL telescope
    )
    # read reflectivities from H5 file
    let resources = parseResourcesPath()
    let h5fname = resources / parseLlnlReflectivityFile()

    let numCoatings = result.layers.len
    var h5f = H5open(h5fname, "r")
    let energies = h5f["/Energy", float]
    let angles = h5f["/Angles", float]
    var reflectivities = newSeq[Interpolator2DType[float]]()
    for i in 0 ..< numCoatings:
      let reflDset = h5f[("Reflectivity" & $i).dset_str]
      let data = reflDset[float].toTensor.reshape(reflDset.shape)
      reflectivities.add newBilinearSpline(
        data,
        (angles.min, angles.max),
        (energies.min, energies.max)
      )
    discard h5f.close()
    result.reflectivities = reflectivities
  of tkXMM:
    result = Reflectivity(
      kind: rkSingleCoating
    )
    # read reflectivities from H5 file
    let resources = parseResourcesPath()
    let h5fname = resources / parseGoldReflectivityFile()

    var h5f = H5open(h5fname, "r")
    let energies = h5f["/Energy", float]
    let angles = h5f["/Angles", float]
    var reflectivities = newSeq[Interpolator2DType[float]]()
    let reflDset = h5f[("Reflectivity").dset_str]
    let data = reflDset[float].toTensor.reshape(reflDset.shape)
    discard h5f.close()

    let spl = newBilinearSpline(
      data,
      (angles.min, angles.max),
      (energies.min, energies.max)
    )
    result.reflectivity = spl
  of tkCustomBabyIAXO, tkAbrixas:
    doAssert false, "Reflectivities are not yet implemented for the " &
      "tkCustomBabyIAXO and tkAbrixas optics."
  else:
    let resources = parseResourcesPath()
    var df = readCsv(resources / parseLlnlTelescopeFile())
    when false: # to be implemented if required. Need the real area of the telescope (e.g. bore
                # diameter, but that's not strictly what it is
      let readRealArea = parseTelescopeArea()
      let inactiveArea = parseInactiveTelArea() # area blocked by opaque structures
      df = df
      .mutate(f{"Transmission" ~ idx("EffectiveArea[cm²]") /
        (readRealArea - inactiveArea)})
    result = Reflectivity(
      kind: rkEffectiveArea,
      telescopeTransmission: newLinear1D(df["Energy[keV]", float].toRawSeq,
                                         df["Transmission", float].toRawSeq)
    )

proc initTelescope(optics: TelescopeKind): Telescope =
  ## TODO: add telescope section to the config file!
  # Note: The required other angles R2, R3, R4 and R5 are computed in the code
  # based on R1 and the angles.
  case optics
  of tkLLNL:
    result = Telescope(
      kind: tkLLNL,
      optics_entrance: @[-83.0, 0.0, 0.0].mapIt(it.mm),
      optics_exit: @[-83.0, 0.0, 454.0].mapIt(it.mm),
      telescope_turned_x: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      telescope_turned_y: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      # the radii of the shells
      allThickness: @[0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                      0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2].mapIt(it.mm),
      # the radii of the shells
      allR1: @[63.006, 65.606, 68.305, 71.105, 74.011, 77.027, 80.157,
               83.405, 86.775, 90.272, 93.902, 97.668, 101.576, 105.632].mapIt(it.mm),
      allXsep: @[4.171, 4.140, 4.221, 4.190, 4.228, 4.245, 4.288, 4.284,
                 4.306, 4.324, 4.373, 4.387, 4.403, 4.481].mapIt(it.mm),
      # the angles of the mirror shells coresponding to the radii above
      allAngles: @[0.579, 0.603, 0.628, 0.654, 0.680, 0.708, 0.737, 0.767,
                   0.798, 0.830, 0.863, 0.898, 0.933, 0.970].mapIt(it.Degree),
      lMirror: 225.0.mm, # Mirror length
      holeInOptics: 0.0.mm, #max 20.9.mm
      numberOfHoles: 5,
      holeType: htCross, #the type or shape of the hole in the middle of the optics
      reflectivity: optics.initReflectivity()
    )
  of tkXMM:
    result = Telescope(
      kind: tkXMM,
      optics_entrance: @[0.0, -0.0, 0.0].mapIt(it.mm),
      optics_exit: @[0.0, -0.0, 600.0].mapIt(it.mm),
      telescope_turned_x: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      telescope_turned_y: 0.0.°, #the angle by which the telescope is turned in respect to the magnet
      # Measurements of the Telescope mirrors in the following, R1 are the radii of the mirror shells at the entrance of the mirror
      allThickness: @[0.468, 0.475, 0.482, 0.490, 0.497, 0.504, 0.511, 0.519, 0.526, 0.534,
                      0.542, 0.549, 0.557, 0.566, 0.574, 0.583, 0.591, 0.600, 0.609, 0.618,
                      0.627, 0.636, 0.646, 0.655, 0.665, 0.675, 0.684, 0.694, 0.704, 0.714,
                      0.724, 0.735, 0.745, 0.756, 0.768, 0.779, 0.790, 0.802, 0.814, 0.826,
                      0.838, 0.850, 0.862, 0.874, 0.887, 0.900, 0.913, 0.927, 0.941, 0.955,
                      0.968, 0.983, 0.997, 1.011, 1.026, 1.041, 1.055, 1.070].mapIt(it.mm),
      # the radii of the shells closest to the magnet, now correct
      allR1: @[153.126, 155.419, 157.731, 160.065, 162.428, 164.812, 167.225, 169.66, 172.124,
               174.608, 177.123, 179.658, 182.224, 184.971, 187.749, 190.556, 193.394, 196.263,
               199.161, 202.09, 205.05, 208.09, 211.16, 214.261, 217.392, 220.553, 223.755,
               226.987, 230.249, 233.552, 236.885, 240.248, 243.652, 247.298, 250.984, 254.711,
               258.478, 262.276, 266.114, 269.992, 273.911, 277.87, 281.869, 285.92, 290.01,
               294.292, 298.676, 303.109, 307.584, 312.108, 316.674, 321.289, 325.955, 330.672,
               335.439, 340.246, 345.104, 350.013].mapIt(it.mm),
      allXsep: @[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].mapIt(it.mm),
      # the angles of the mirror shells coresponding to the radii above, now correct
      allAngles: @[0.29, 0.294, 0.298, 0.303, 0.307, 0.312, 0.316, 0.321, 0.325, 0.33, 0.335,
                   0.34, 0.345, 0.35, 0.355, 0.36, 0.366, 0.371, 0.377, 0.382, 0.388, 0.393,
                   0.399, 0.405, 0.411, 0.417, 0.423, 0.429, 0.435, 0.441, 0.448, 0.454, 0.461,
                   0.467, 0.474, 0.481, 0.489, 0.496, 0.503, 0.51, 0.518, 0.525, 0.533, 0.54,
                   0.548, 0.556, 0.564, 0.573, 0.581, 0.59, 0.598, 0.607, 0.616, 0.625, 0.634,
                   0.643, 0.652, 0.661].mapIt(it.Degree),
      lMirror: 300.0.mm, # Mirror length
      holeInOptics: 0.2.mm, #max 20.9.mm
      numberOfHoles: 1,
      holeType: htNone, #the type or shape of the hole in the middle of the optics
      reflectivity: optics.initReflectivity()
    )
  else:
    doAssert false, "The telescope for kind " & $optics & " has not been implemented yet!"

proc initTestXraySource(setup: ExperimentSetupKind, flags: set[ConfigFlags]): TestXraySource =
  let active = cfXrayTest in flags
  let sourceOpt = maybeParseTestXraySource(flags)
  if sourceOpt.isSome:
    return sourceOpt.get
  case setup
  of esCAST:
    result = TestXraySource(
      active: active,
      distance: 100.0.mm, #distance between the entrance of the magnet an a test Xray source
      radius: 10.0.mm,
      offAxisUp: 200.0.mm,
      offAxisLeft: 0.0.mm,
      lengthCol: 50.0.mm,
      energy: 1.0.keV,
      activity: 1.0.GBq
    )
  of esBabyIAXO:
    result = TestXraySource(
      active: active,
      distance: 2000.0.mm, #88700.0.mm, #distance between the entrance of the magnet an a test Xray source
      radius: 350.0.mm,
      offAxisUp: 0.0.mm,
      offAxisLeft: 0.0.mm,#4.5.mm, #
      lengthCol: 0.0.mm,
      energy: 0.021.keV,
      activity: 0.125.GBq, #proposed source thing by Thomas #1.0.GBq,
    )

proc initDetectorInstallation(setup: ExperimentSetupKind,
                              flags: set[ConfigFlags]): DetectorInstallation =
  let detOpt = maybeParseDetectorInstallation(flags)
  if detOpt.isSome:
    return detOpt.get
  case setup
  of esCAST:
    result = DetectorInstallation(
      distanceDetectorXRT: 1485.0.mm, #1300.0 # is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 # was the Telescope of 2014 (MPE XRT) also: Aperatur changed #ok
      distanceWindowFocalPlane: 0.0.mm, # #no change, because don't know

      lateralShift: 0.0.mm, #lateral ofset of the detector in repect to the beamline
      transversalShift: 0.0.mm #transversal ofset of the detector in repect to the beamline #0.0.mm #
    )
  of esBabyIAXO:
    result = DetectorInstallation(
      distanceDetectorXRT: 7500.0.mm, # # one possibility, the other is 5050 mm
      distanceWindowFocalPlane: 0.0.mm, # #no change, because don't know #good idea
      lateralShift: 0.0.mm, #(sin(0.0.degToRad) * 7500.0).mm, #lateral ofset of the detector in repect to the beamline #0.0.mm #
      transversalShift: (sin(0.0.degToRad) * 7500.0).mm #-0.0.mm # ##transversal ofset of the detector in repect to the beamline #0.0.mm #
    )


proc getOpticForSetup(setup: ExperimentSetupKind): TelescopeKind =
  ## XXX: replace this by user input & config file selection. Only fall back if
  ## neither provided!
  case setup
  of esCAST: result = tkLLNL
  of esBabyIAXO: result = tkXMM # tkCustomBabyIAXO

proc newExperimentSetup*(setup: ExperimentSetupKind,
                         stage: StageKind,
                         # optics: TelescopeKind,
                         flags: set[ConfigFlags]): ExperimentSetup =
  let optics = getOpticForSetup(setup)
  result = ExperimentSetup(
    kind: setup,
    stage: stage,
    magnet: setup.initMagnet(flags),
    telescope: optics.initTelescope(),
    testSource: setup.initTestXraySource(flags),
    pipes: setup.initPipes(),
    detectorInstall: setup.initDetectorInstallation(flags)
  )

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
  result.detectorWindowAperture = ChipXMax
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

proc computeReflectivity(expSetup: ExperimentSetup, energy: keV,
                         hitLayer: int, # the layer of the telescope hit
                         transmissionMagnet: float,
                         p, ya: float, alpha1, alpha2: Degree,
                         flags: set[ConfigFlags]): tuple[reflect: float,
                                                         weight: float] =
  ## Computes the reflectivity of the two reflections in the telescope of `expSetup`
  ## according to the given `angles1/2` and the `energy`.
  if cfIgnoreReflection in flags: # short circut for
    result.reflect = 1.0
    result.weight = transmissionMagnet
    return result

  let refl = expSetup.telescope.reflectivity
  case refl.kind
  of rkEffectiveArea:
    # we compute the pitch and yaw angles under which the ray reaches the telescope.
    # Based on a very coarse set of data points for the efficiency under yaw/pitch angles
    # a polynomial approximation (for the LLNL CAST telescope!) is used to compute the
    # final efficiency.
    let
      transmissionTelescopePitch = (0.0008*p*p*p*p + 1e-04*p*p*p - 0.4489*p*p -
          0.3116*p + 96.787) / 100.0
      transmissionTelescopeYaw = (6.0e-7 * pow(ya, 6.0) - 1.0e-5 * pow(ya, 5.0) -
          0.0001 * pow(ya, 4.0) + 0.0034 * pow(ya, 3.0) - 0.0292 * pow(ya, 2.0) -
          0.1534 * ya + 99.959) / 100.0
    let transmissionTel = refl.telescopeTransmission.eval(energy)
    result.reflect = transmissionTel * transmissionTelescopePitch * transmissionTelescopeYaw
    # transmission probabilities times axion emission rate times the flux fraction
    result.weight = (result.reflect * transmissionMagnet)
  of rkSingleCoating:
    ## TODO: This needs to be replaced by a 2D interpolation!
    {.gcsafe.}: # `{.gcsafe.}` due to closures in `numericalnim` missing the tag
      let
        reflectionProb1 = refl.reflectivity.eval(alpha1.float, energy.float)
        reflectionProb2 = refl.reflectivity.eval(alpha2.float, energy.float)
      result.reflect = reflectionProb1 * reflectionProb2
      result.weight = result.reflect * transmissionMagnet
  of rkMultiCoating:
    # use the hit layer to know which interpolator we have to use
    let layerIdx = refl.layers.lowerBound(hitLayer)
    let reflLayer = refl.reflectivities[layerIdx]
    {.gcsafe.}: # `{.gcsafe.}` due to closures in `numericalnim` missing the tag
      let
        reflectionProb1 = reflLayer.eval(alpha1.float, energy.float)
        reflectionProb2 = reflLayer.eval(alpha2.float, energy.float)
      result.reflect = reflectionProb1 * reflectionProb2
      result.weight = result.reflect * transmissionMagnet

proc computeMagnetTransmission(
  expSetup: ExperimentSetup, energy: keV, distancePipe: Meter, pathCB: mm, ya: float,
  flags: set[ConfigFlags]
     ): float =
  ## Computes the effective 'transmission' through the magnet. This means the conversion
  ## probability for axions to X-rays and in case of a buffer gas in the cold bore also the
  ## intensity suppression of the generated X-rays after conversion.
  ##
  ## TODO: the combination of conversion probability and intensity suppression of X-rays
  ## in a gas, is more complicated than this!
  case expSetup.stage
  of skVacuum:
    let probConversionMagnet = if cfIgnoreConvProb notin flags:
                                 conversionProb(expSetup.magnet.B, g_aγ, pathCB)
                               else:
                                 1.0
    result = cos(ya) * probConversionMagnet
  of skGas:
    let
      pGas = expSetup.magnet.pGasRoom / RoomTemp * expSetup.magnet.tGas
      ## TODO: use `unchained` in `axionMass` and avoid manual float conversions / use `to`
      effPhotonMass = effPhotonMass2(
        pGas.float, pathCB.to(m).float,
        expSetup.magnet.radiusCB.to(m).float,
        expSetup.magnet.tGas.float
      )
      probConversionMagnetGas =
        if cfIgnoreConvProb notin flags:
          axionConversionProb2(
            mAxion, energy, pGas.float,
            expSetup.magnet.tGas.float, pathCB.to(m).float,
            expSetup.magnet.radiusCB.to(m).float,
            g_aγ.float, expSetup.magnet.B.float
          ) # for setup including gas: functions are in axionmass/axionMassforMagnet
        else:
          1.0 # else perfect transmission
      absorbtionXrays = intensitySuppression2(
        energy, pathCB.to(m).float,
        distancePipe.float,
        pGas.float,
        expSetup.magnet.tGas.float,
        RoomTemp.float) #room temperature in K
    #echo "axion mass in [eV] ", mAxion, " effective photon mass in [eV] ", effPhotonMass
    result = cos(ya) * probConversionMagnetGas * absorbtionXrays # for setup with gas

func radiusAndPhi(vec: Vec3[float]): (MilliMeter, float) =
  ## Simple helper that computes the radius and angle φ from the given vector
  ## in the x/y plane
  let
    radius = sqrt((vec[0].mm) * (vec[0].mm) + (vec[1].mm) * (vec[1].mm))
    phi = radtoDeg(arccos(vec[0].mm / radius))
  result = (radius, phi)

proc lineIntersectsOpaqueTelescopeStructures(
  expSetup: ExperimentSetup,
  radialDist: MilliMeter, testVector, vectorXRT, pointExitCB, pointEntranceXRT: Vec3[float],
  onlyInnerPart = false
     ): bool =
  ## Checks whether the incoming ray hits any of the fully opaque parts of the telescope,
  ## usually the support structure holding the shells in place.
  case expSetup.telescope.kind
  of tkLLNL:
    ## there is a 2mm wide graphite block between each glass mirror, to seperate them
    ## in the middle of the X-ray telescope. Return if hit
    if pointEntranceXRT[1] <= 1.0 and pointEntranceXRT[1] >= -1.0: return
  of tkXMM:
    ## BabyIAXO return if X-ray its the spider structure
    ## here we have a spider structure for the XMM telescope:
    ### 3D spider structure 1st draft ###

    var factorSpider = (-85.0 - pointExitCB[2]) / vectorXRT[2]
    var pointEntranceSpider = pointExitCB + factorSpider * vectorXRT
    let
      (radius, phiFlat) = radiusAndPhi(pointEntranceXRT)
      (radiusSpider, phiFlatSpider) = radiusAndPhi(pointEntranceSpider)
    # `onlyInnerPart` is used later in code to force entering this branch
    if onlyInnerPart or radialDist <= 64.7.mm:
      let nHoles = expSetup.telescope.numberOfHoles
      for l in -(nHoles - ceil(nHoles.float / 2.0).int) .. (nHoles - ceil(nHoles.float / 2.0).int):
        var centerHole = testVector
        if l != 0:
          if abs(l) mod 2 == 0:
            centerHole[1] += 2.0 * l.float * expSetup.telescope.holeInOptics.float
          else:
            centerHole[0] += 2.0 * (l + (l / abs(l))).float * expSetup.telescope.holeInOptics.float
        if lineIntersectsObject(expSetup.telescope.holeType, pointExitCB, pointEntranceXRT,
                                centerHole, expSetup.telescope.holeInOptics):
          result = false
          break
        else:
          result = true
    # hole max 20.9 mm doesnt really matter because there are no mirrors in the
    # middle and these axions dont reach the window anyways
    elif radialDist < 151.6.mm and radialDist > (151.6 - 20.9).mm:
      result = true
    elif radialDist > 64.7.mm:
      # iterate all strips of the spider structure
      for i in 0 .. 16:
        # spider strips (actually wider for innermost but doesn't matter because
        # it doesnt reach the window anyways)
        if ((phiFlat >= (-1.145 + 22.5 * i.float) and phiFlat <= (1.145 + 22.5 * i.float))) or
           ((phiFlatSpider >= (-1.145 + 22.5 * i.float) and phiFlatSpider <= (1.145 + 22.5 * i.float))):
          result = true
          break
  else:
    doAssert false, "The telescope kind " & $expSetup.telescope.kind & " does not have any opaque " &
      "structures implemented yet."

proc lineHitsNickel(expSetup: ExperimentSetup, α1: Degree, r1: MilliMeter,
                    hitLayer: int, # the shell of the telescope hit
                    pointMirror1, pointAfterMirror1, pointMirror2: Vec3[float]): bool =
  ## Checks if the ray described given `α1` (the angle from the optical axis)
  ## hits the nickel after reflection
  ##
  ## Note: to generalize this to the "second part" of the telescope the arguments in
  ## principle remain almost unchanged (`α1` becomes a generic `α` which gets `angle2`
  ## from the calling code and the code uses `pointMirror2` instead of `pointMirror1`.
  ## Further the radius `r1` must become the equivalent `r4` which would have to be computed
  ## and the lengths adjusted.
  const verbose = false
  template tel(): untyped = expSetup.telescope
  if hitLayer > 0: # cannot hit nickel if ray hits lowest shell (no nickel below!)
    let tanα = tan(α1.to(Radian))
    let hL = hitLayer - 1 # layer below the hit layer
    let compVal = (r1 - (tel.allR1[hL] + tel.allThickness[hL])) /
       (tel.lMirror - pointMirror1[2].mm)
    when verbose:
      let pointLowerMirror = findPosXRT(
        pointAfterMirror1, pointMirror1,
        tel.allR1[hL] + tel.allThickness[hL],
        tel.allR1[hL] -
          tel.lMirror * sin(tel.allAngles[hL].to(Radian)) +
          tel.allThickness[hL],
        tel.allAngles[hL].to(Radian),
        tel.lMirror,
        0.0.mm, 0.01.mm, 0.0.mm, 2.5.mm
      )
      echo "hit Nickel"
      echo pointLowerMirror, " ", pointMirror1, " ", pointMirror2, " ", angle1[1], " ", angle2[1], " ", h
    result = tanα > compVal

proc traceAxion(res: var Axion,
                centerVecs: CenterVectors,
                expSetup: ExperimentSetup,
                detectorSetup: DetectorSetup,
                fluxRadiusCDF: seq[float],
                diffFluxCDFs: seq[seq[float]],
                energies: seq[keV],
                flags: set[ConfigFlags]
               ) =
  ## Check if we run with an X-ray source or compute from the Sun
  let testXray = expSetup.testSource.active

  var rayOrigin: Vec3[float] # Starting point of the ray (typically axion somewhere in the Sun)
  var pointExitCBMagneticField: Vec3[float] # position at exit of magnet (entry from viewpoint of incoming ray)
  var energyAx: keV # energy of the axion
  if not testXray:
    # Get a random point in the sun, biased by the emission rate, which is higher
    # at smalller radii, so this will give more points in the center of the sun
    rayOrigin = getRandomPointFromSolarModel(centerVecs.sun, RadiusSun, fluxRadiusCDF)
    # Get a random point at the end of the coldbore of the magnet to take all
    # axions into account that make it to this point no matter where they enter the magnet
    pointExitCBMagneticField = getRandomPointOnDisk(
      centerVecs.exitCBMagneticField,
      expSetup.magnet.radiusCB
    )
    # Get a random energy for the axion biased by the emission rate ##
    energyAx = getRandomEnergyFromSolarModel(
      rayOrigin, centerVecs.sun, RadiusSun, energies, diffFluxCDFs
    )
  else:
    ## XXX: CLEAN THIS UP! likely just remove dead code & create procs for the "complicated parts" that
    ## yield a single value for `pointExitCBMagneticField`!
    rayOrigin = getRandomPointOnDisk(
      centerVecs.xraySource, expSetup.testSource.radius
    )
    energyAx = expSetup.testSource.energy

    ## i.e. all this here
    var centerSpot = vec3(0.0)
    if abs(expSetup.testSource.offAxisUp.float) > 50.0:
      centerSpot[0] = expSetup.testSource.offAxisLeft.float
      centerSpot[1] = expSetup.testSource.offAxisUp.float
    centerSpot[2] = expSetup.magnet.lengthB.float
    #var radiusProjection = expSetup.testSource.radius * 2.0 * (- centerVecs.xraySource[2] - expSetup.lengthCol.float).mm / expSetup.lengthCol + expSetup.radiusXraySource
    var radiusSpot = expSetup.telescope.holeInOptics.float
    case expSetup.telescope.holeType
    of htCross, htStar:
      radiusSpot *= (expSetup.telescope.numberOfHoles.float + 16.0)
    of htNone:
      radiusSpot = expSetup.testSource.radius.float / 100.0 #* 3.0 #expSetup.telescope.radiusCB.float / 6.0 #
    else:
      ## XXX: check if this branch should be triggered sometimes. What hole type does this correspond to?
      radiusSpot *= (expSetup.telescope.numberOfHoles.float + 3.0)
    pointExitCBMagneticField[0] = rayOrigin[0] + rand(0.1) - 0.05 #for parallel light
    pointExitCBMagneticField[1] = rayOrigin[1] + rand(0.1) - 0.05 #for parallel light
    ## XXX: this is my modification, as this is the value that would be at [2] if we still sampled from
    ## via `getRandomPointOnDisk` before entering this `else` branch. CHECK THIS!
    pointExitCBMagneticField[2] = expSetup.magnet.lengthB.float

    #pointExitCBMagneticField = getRandomPointOnDisk(centerSpot, (radiusSpot).mm) # for more statistics with hole through optics
    if not lineIntersectsCircle(rayOrigin, pointExitCBMagneticField, centerVecs.collimator, expSetup.testSource.radius):
      return

    var xraysThroughHole = PI * radiusSpot * radiusSpot /
      (4.0 * PI * (- centerVecs.xraySource[2] + centerVecs.exitPipeVT3XRT[2]).mm *
      (- centerVecs.xraySource[2] + centerVecs.exitPipeVT3XRT[2]).mm) * expSetup.testSource.activity
    var testTime = 1_000_000 / (xraysThroughHole * 3600.0.s * 24.0)
    #echo "Days Testing ", testTime, " with a ", expSetup.testSource.activity, " source"
  #let emissionRateAx = getRandomEmissionRateFromSolarModel(
  #  rayOrigin, centerVecs.centerSun, RadiusSun, emRates, emRateCDFs
  #)
  ## Throw away all the axions, that don't make it through the piping system and therefore exit the system at some point ##


  let intersectsEntranceCB = lineIntersectsCircle(rayOrigin,
      pointExitCBMagneticField, centerVecs.entranceCB, expSetup.magnet.radiusCB)
  var intersectsCB = false
  var rayOriginInSun = rayOrigin - centerVecs.sun #
  #echo energyAx.float, " ", sqrt(rayOriginInSun[0].float * rayOriginInSun[0].float + rayOriginInSun[1].float * rayOriginInSun[1].float + rayOriginInSun[2].float * rayOriginInSun[2].float) / 6.9e11, " ", energyAx.float  * energyAx.float * (rayOriginInSun[0].float * rayOriginInSun[0].float + rayOriginInSun[1].float * rayOriginInSun[1].float + rayOriginInSun[2].float * rayOriginInSun[2].float) / 6.9e11 / 6.9e11
  res.emratesPre = 1.0 #* energyAx.float * energyAx.float * (rayOriginInSun[0].float * rayOriginInSun[0].float + rayOriginInSun[1].float * rayOriginInSun[1].float + rayOriginInSun[2].float * rayOriginInSun[2].float) / 6.9e11 / 6.9e11
  res.energiesPre = energyAx
  if (not intersectsEntranceCB):

    intersectsCB = lineIntersectsCylinderOnce(rayOrigin,
        pointExitCBMagneticField, centerVecs.entranceCB,
        centerVecs.exitCB, expSetup.magnet.radiusCB)
  if (not intersectsEntranceCB and not intersectsCB): return

  var intersect = vec3(0.0) #isnt't changed for axions that hit the entrance of the coldbore because the z value is 0 then anyways
  if (not intersectsEntranceCB): #generates problems with the weight because the weight is multiplied with the difference of the leght of the path of the particle and the legth of the coldbore
    intersect = getIntersectLineIntersectsCylinderOnce(
      rayOrigin,
      pointExitCBMagneticField, centerVecs.entranceCB,
      centerVecs.exitCB, expSetup.magnet.radiusCB
    ) #rayOrigin + ((centerVecs.entranceCB[2] - rayOrigin[2]) / (pointExitCBMagneticField - rayOrigin)[2]) * (pointExitCBMagneticField - rayOrigin)
    #echo "doesn't intersect magnet entrance: ", intersect
  else:
    intersect = getIntersectlineIntersectsCircle(
      rayOrigin,
      pointExitCBMagneticField, centerVecs.entranceCB
    )
    #echo "does intersect magnet entrance: ", intersect

  ##get the length of the path of the axion in the magnetic field to get the probability of conversion later
  let pathCB = (pointExitCBMagneticField - intersect).length.mm

  ## Return if... WRITE ME
  if not lineIntersectsCircle(rayOrigin, pointExitCBMagneticField,
                              centerVecs.exitCB, expSetup.magnet.radiusCB):
    return

  var pointExitCB = rayOrigin +
    ((centerVecs.exitCB[2] - rayOrigin[2]) /
     (pointExitCBMagneticField - rayOrigin)[2]) *
    (pointExitCBMagneticField - rayOrigin)

  ## Return if... WRITE ME
  if not lineIntersectsCircle(pointExitCBMagneticField, pointExitCB,
                              centerVecs.exitPipeCBVT3, expSetup.pipes.coldBoreToVT3.radius):
    return

  let pointExitPipeCBVT3 = pointExitCBMagneticField +
   ((centerVecs.exitPipeCBVT3[2] - pointExitCBMagneticField[2]) /
    (pointExitCB - pointExitCBMagneticField)[2]) *
   (pointExitCB - pointExitCBMagneticField)

  ## Return if... WRITE ME
  if not lineIntersectsCircle(pointExitCB, pointExitPipeCBVT3,
                              centerVecs.exitPipeVT3XRT, expSetup.pipes.coldBoreToVT3.radius):
    return

  let pointExitPipeVT3XRT = pointExitCB +
    ((centerVecs.exitPipeVT3XRT[2] - pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) *
    (pointExitPipeCBVT3 - pointExitCB)

  var vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB

  #echo centerVecs.collimator, " ", pointExitPipeVT3XRT
  ###################from the CB (coldbore(pipe in Magnet)) to the XRT (XrayTelescope)#######################
  var vectorXRT = vectorBeforeXRT
  let
    turnedX = expSetup.telescope.telescope_turned_x.to(Radian)
    turnedY = expSetup.telescope.telescope_turned_y.to(Radian)
  vectorXRT = rotateInY(rotateInX(vectorXRT, turnedX), turnedY)

  pointExitCB[2] -= centerVecs.exitPipeVT3XRT[2]
  pointExitCB = rotateInY(rotateInX(pointExitCB, turnedX), turnedY) -
                vec3(expSetup.telescope.optics_entrance[0].float, expSetup.telescope.optics_entrance[1].float, 0.0)

  let
    factor = (0.0 - pointExitCB[2]) / vectorXRT[2]
    pointEntranceXRT = pointExitCB + factor * vectorXRT
  vectorBeforeXRT = vectorXRT
  ## Coordinate transform from cartesian to polar at the XRT entrance
  let
    # compute the radial distance of the ray at the entrance of the XRT to the optical
    # axis, as this is the basis to determine which layer will be hit (by comparing with
    # the height of the layer)
    (radialDist, _) = radiusAndPhi(pointEntranceXRT)

  ## Check if any of the opaque structures of the telescopes are hit (graphite block for LLNL,
  ## spider structure for XMM, ...)
  if expSetup.lineIntersectsOpaqueTelescopeStructures(
    radialDist, vec3(0.0), vectorXRT, pointExitCB, pointEntranceXRT
  ):
    return

  ## Calculate the way of the axion through the telescope by manually reflecting the ray on the two mirror layers and then ending up before the detector ##

  var
    minDist: MilliMeter = Inf.mm # the minimal (positive) distance between ray and telescope layers
    r1 = 0.0.mm
    r2 = 0.0.mm
    r3 = 0.0.mm
    r4 = 0.0.mm
    r5 = 0.0.mm
    beta = 0.0.Radian ## in degree
    xSep = 0.0.mm
    hitLayer: int # the integer of the layer that was hit
    lengthTelescope = (expSetup.telescope.lMirror + 0.5 * xSep) * cos(expSetup.telescope.allAngles[0].to(Radian)) +
                      (expSetup.telescope.lMirror + 0.5 * xSep) * cos(3.0 * expSetup.telescope.allAngles[0].to(Radian))
    centerEndXRT = vec3(0.0)
  centerEndXRT[0] = 0.0 #lengthTelescope.float * tan(expSetup.telescope_turned.to(Radian)) thats not corrrect but the thing is already turned
  centerEndXRT[2] = lengthTelescope.float

  let allR1 = expSetup.telescope.allR1
  # remove rays that are further out than outer most layer
  if radialDist > allR1[allR1.high]: return
  # walk all layers and check which one was hit, if glass is hit remove ray
  ## TODO: possible optimization: as the layers are ordered from inner to outer, don't
  ## we always go from most negative to most positive in distance? I.e. the first positive
  ## distance value is also the smallest one? If so, we could immediately break after that
  ## layer!
  for j in 0 ..< allR1.len:
    # get rid of where the X-rays hit the glass frontal
    if radialDist > allR1[j] and
        radialDist < (allR1[j] + expSetup.telescope.allThickness[j]):
      return
    let dist = allR1[j] - radialDist # distance of ray from layer j
    # if we are not above the layer and smaller than latest `minDist`, assign
    if dist > 0.0.mm and dist < minDist:
      minDist = dist
      hitLayer = j # j-th layer currently hit layer
      # now get/compute the radii for the hit layer
      r1 = allR1[j]
      beta = expSetup.telescope.allAngles[j].to(Radian)
      xSep = expSetup.telescope.allXsep[j]
      r2 = r1 - expSetup.telescope.lMirror * sin(beta)
      r3 = r2 - 0.5 * xSep * tan(beta)
      r4 = r3 - 0.5 * xSep * tan(3.0 * beta)
      r5 = r4 - expSetup.telescope.lMirror * sin(3.0 * beta)

      ## TODO: verify that `break` here is actually justified
      # break

  ## Why do we check for `lineIntersectsOpaqueTelescopeStructure` here again?
  if testXray and minDist > 100.0.mm and
     expSetup.lineIntersectsOpaqueTelescopeStructures(
       radialDist, centerEndXRT, vectorXRT, pointExitCB, pointEntranceXRT,
       onlyInnerPart = true
     ):
    return

  let
    beta3 = 3.0 * beta
    distanceMirrors = cos(beta) * (xSep + expSetup.telescope.lMirror)
    pointMirror1 = findPosXRT(pointEntranceXRT, pointExitCB, r1,
                              r2, beta, expSetup.telescope.lMirror, 0.0.mm, 0.001.mm, 1.0.mm, 1.1.mm, "Mirror 1")


  let
    vectorAfterMirror1 = getVectoraAfterMirror(pointEntranceXRT,
        pointExitCB, pointMirror1, beta, "vectorAfter")
    pointAfterMirror1 = pointMirror1 + 200.0 * vectorAfterMirror1
    pointMirror2 = findPosXRT(
      pointAfterMirror1, pointMirror1, r4, r5, beta3,
      expSetup.telescope.lMirror, distanceMirrors, 0.01.mm, 0.0.mm, 2.5.mm, "Mirror 2"
    )


  #echo r1, " ", allR1[hitLayer - 1], " ", r2, " ", allR1[hitLayer - 1] - expSetup.telescope.lMirror * sin(expSetup.allAngles[hitLayer - 1].to(Radian))


  if pointMirror2[0] == 0.0 and pointMirror2[1] == 0.0 and pointMirror2[2] ==
      0.0: return ## with more uncertainty, 10% of the 0.1% we loose here can be recovered, but it gets more uncertain
  let
    vectorAfterMirrors = getVectoraAfterMirror(
      pointAfterMirror1, pointMirror1, pointMirror2, beta3, "vectorAfter"
    )
    pointAfterMirror2 = pointMirror2 + 200.0 * vectorAfterMirrors
  let
    angle1 = getVectoraAfterMirror(
      pointEntranceXRT,
      pointExitCB, pointMirror1, beta, "angle"
    )
    angle2 = getVectoraAfterMirror(
      pointAfterMirror1, pointMirror1, pointMirror2, beta3, "angle"
    )
    # `getVectoraAfterMirror` with `"angle"` argument returns `Degree`
    alpha1 = angle1[1].Degree
    alpha2 = angle2[1].Degree

  # getting rid of the X-rays that hit the shell below
  res.hitNickel = expSetup.lineHitsNickel(
    alpha1, r1, hitLayer,
    pointMirror1, pointAfterMirror1, pointMirror2
  )
  # if nickel was hit, stop here
  if res.hitNickel:
    return

  ## Check if `findPosXRT` did not find a hit on the layer. In that case the two `pointMirror`
  ## variables actually contain the same (at least `z`) data, because the input is returned
  ## unchanged in that case.
  let z0 = pointExitCB[2]
  let z1 = pointMirror1[2]
  let z2 = pointMirror2[2]
  if almostEqual(z1, z2) or almostEqual(z1, z0):
    return
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
    distDet = distanceMirrors - 0.5 * expSetup.telescope.allXsep[8] * cos(beta) +
        expSetup.detectorInstall.distanceDetectorXRT -
        expSetup.detectorInstall.distanceWindowFocalPlane # distance of the detector from the entrance of the optics
    d = - expSetup.telescope.optics_entrance[0]
  pointDetectorWindow = getPointDetectorWindow(
    pointMirror2,
    pointAfterMirror2,
    distDet, expSetup.telescope.allXsep[8], d, expSetup.pipes.pipesTurned
  )
  pointEndDetector = getPointDetectorWindow(
    pointMirror2, pointAfterMirror2,
    (distDet + detectorSetup.depthDet),
    expSetup.telescope.allXsep[8], d, expSetup.pipes.pipesTurned
  )

  res.deviationDet = sqrt(
    pow((pointEndDetector[0] - pointDetectorWindow[0]), 2.0) +
    pow((pointEndDetector[1] - pointDetectorWindow[1]), 2.0)
  )

  var valuesPix = getPixelValue(pointEntranceXRT)
  res.pointdataXBefore = pointEntranceXRT[0]
  res.pointdataYBefore = pointEntranceXRT[1]
  res.pixvalsX = valuesPix[0]
  res.pixvalsY = valuesPix[1]

  ## Calculate the weight for each axion as the probability of its signal arriving at the detector which consists of:
  ## The probability of conversion into Xrays in the path in the homogeneous magnetic field
  ## NOT! The fraction of the flux of actions emitted from the given point in the sun per second, that actually arrives at the random point in the coldbore end, because this is handled by biasing the random origin
  ## The transmission probability through the LLNL Xray telescope, which depends on the energy of the Xray and its angle

  vectorBeforeXRT = - vectorBeforeXRT # because we have to get the angles from the perspective of the XRT

  let vecLength = vectorBeforeXRT.length
  var vectorBeforeXRTPolar = vec3(
    vecLength, # r
    radToDeg(arccos(vectorBeforeXRT[0] / vecLength)), # θ
    radToDeg(arctan2(vectorBeforeXRT[2], vectorBeforeXRT[1])) # φ
  )

  #this is the pitch angle # not sure why plus 90
  vectorBeforeXRTPolar[1] = vectorBeforeXRTPolar[1] - 90.0
  let p = vectorBeforeXRTPolar[1]
  #this is the yaw angle, floor to roof
  vectorBeforeXRTPolar[2] = vectorBeforeXRTPolar[2] + 90.0
  let ya = vectorBeforeXRTPolar[2]
  let distancePipe = (pointDetectorWindow[2] - pointExitCB[2]).mm.to(m) # needs to be in meter

  ## this is the transformation probability of an axion into a photon, if an axion
  ## flying straight through the magnet had one of 100%, angular dependency of the primakoff effect
  res.transmissionMagnet = expSetup.computeMagnetTransmission(
    energyAx, distancePipe, pathCB, ya, flags
  )
  res.yawAngles = ya

  var weight = 1.0
  (res.reflect, weight) = expSetup.computeReflectivity(
    energyAx, hitLayer, res.transmissionMagnet, p, ya, alpha1, alpha2, flags
  )

  if testXray != false:
    if cfXrayTest notin flags:
      weight = res.transmissionMagnet
    n = (distDet - pointExitCB[2].mm) / (pointEntranceXRT - pointExitCB)[2].mm
    pointDetectorWindow = pointExitCB + n.float * (pointEntranceXRT - pointExitCB)
  pointDetectorWindow[0] -= expSetup.detectorInstall.lateralShift.float
  pointDetectorWindow[1] -= expSetup.detectorInstall.transversalShift.float
  if weight != 0:
    res.passedTillWindow = true

  ##Detector window:##
  if cfIgnoreDetWindow notin flags and sqrt(
      pointDetectorWindow[0].mm * pointDetectorWindow[0].mm +
      pointDetectorWindow[1].mm * pointDetectorWindow[1].mm
      ) > detectorSetup.radiusWindow:
      return

  else:
    if abs(pointDetectorWindow[0].mm) > ChipCenterX or abs(pointDetectorWindow[1].mm) > ChipCenterY:
      return

  let
    theta = detectorSetup.theta
    pointDetectorWindowTurned = rotateAroundZ(pointDetectorWindow, theta)
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
  ## XXX: move to a separate proc!
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
  res.shellNumber = hitLayer

  ###detector COS has (0/0) at the bottom left corner of the chip
  pointRadialComponent = sqrt(pointDetectorWindow[0]*pointDetectorWindow[0]+pointDetectorWindow[1]*pointDetectorWindow[1])
  res.pointdataR = pointRadialComponent
  pointDetectorWindow[0] = - pointDetectorWindow[0] + ChipCenterX.float # for the view from the detector to the sun
  pointDetectorWindow[1] = pointDetectorWindow[1] + ChipCenterY.float

  ## TODO: replace by calculation based on user input solar tracking time
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

  ## finally set the `passed` field to indicate this axion went all the way
  if weight != 0:
    res.passed = true

proc traceAxionWrapper(axBuf: ptr UncheckedArray[Axion],
                       bufLen: int,
                       centerVecs: CenterVectors,
                       expSetup: ExperimentSetup,
                       detectorSetup: DetectorSetup,
                       fluxRadiusCDF: seq[float],
                       diffFluxCDFs: seq[seq[float]],
                       energies: seq[keV],
                       flags: set[ConfigFlags]
                       ) =
  echo "Starting weave!"
  parallelFor iSun in 0 ..< bufLen:
    captures: { axBuf, centerVecs, expSetup, fluxRadiusCDF, diffFluxCDFs,
                detectorSetup,
                energies,
                flags }
    axBuf[iSun].traceAxion(centerVecs,
                           expSetup,
                           detectorSetup,
                           fluxRadiusCDF, diffFluxCDFs,
                           energies,
                           flags)

proc generateResultPlots(axions: seq[Axion],
                         windowYear: WindowYearKind,
                         suffix = ""
                         ) =
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
  extractPass(transmissionMagnet)
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
  if dfRad.len > 0:
    ggplot(dfRad, aes("Radial component [mm]", weight = "Transmission probability")) +
      geom_histogram(binWidth = 0.001) +
      ggtitle("Radial distribution of the axions") +
      ggsave(&"../out/radialDistribution_{windowYear}.pdf")

  when false:
    let dfFluxE = seqsToDf({ "Axion energy [keV]": energiesAx.mapIt(it.float),
                             "Transmission probability": weights })

    ggplot(dfFluxE, aes("Axion energy [keV]", weight = "Transmission probability")) +
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
      ggsave(&"../out/fluxAfter_{windowYear}.png")

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

  when false:
    let dfMag = seqsToDf({"Transmission probability": transmissionMagnet,
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
      ggsave(&"../out/transMagnetE_{windowYear}.pdf")

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
  let dfRadFilter = dfRadSig.filter(f{`x` < (ChipCenterX.float + 0.05)}).filter(f{`x` > (ChipCenterX.float - 0.05)})
  if dfRadFilter.len > 0:
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
  dfFluxE2.write_csv(&"axion_gae_1e13_gagamma_{g_aγ.float}_flux_after_exp_N_{NumberOfPointsSun}.csv")
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
    endX = ChipXMax.float  #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endY = ChipYMax.float  #- distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareHeatmap(3000, 3000, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights,
      NumberOfPointsSun.float) #colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareHeatmap(256, 256, beginX, endX, beginY, endY,
      pointdataX, pointdataY, weights, 1.0)
  var heatmaptable3 = prepareHeatmap(3000, 3000, beginX, endX, beginY, endY,
                                     pointdataX, pointdataY, weights, heatmaptable2.max) # if change number of rows: has to be in the maxVal as well
 # echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
 # echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]
  plotHeatmap("Axion Model Fluxfraction", heatmaptable2, 256, $windowYear, rSigma1W, rSigma2W, suffix) #rSigma1, rSigma2)

proc initFullSetup(setup: ExperimentSetupKind,
                   detectorSetup: DetectorSetupKind,
                   stage: StageKind,
                   flags: set[ConfigFlags]): FullRaytraceSetup =
  let expSetup = newExperimentSetup(setup, stage, flags)

  ## TODO: make the code use tensor for the emission rates!
  let resources = parseResourcesPath()
  var emRatesDf = readCsv(resources / parseSolarModelFile())

  # get all radii and energies from DF so that we don't need to compute them manually (risking to
  # messing something up!)
  # sort both just to make sure they really *are* in ascending order
  let radii = emRatesDf["Radius"]
    .unique()
    .toTensor(float)
    .toSeq1D
    .sorted(SortOrder.Ascending)
  let energies = emRatesDf["Energy [keV]"]
    .unique()
    .toTensor(float)
    .toSeq1D
    .mapIt(it.keV)
    .sorted(SortOrder.Ascending)
  var emRates = newSeq[seq[float]]()
  ## group the "solar model" DF by the radius & append the emission rates for all energies
  ## to the `emRates`
  for tup, subDf in groups(emRatesDf.group_by("Radius")):
    doAssert subDf["Energy [keV]", float].toSeq1D.mapIt(it.keV) == energies
    emRates.add subDf["emRates", float].toSeq1D

  var
    fluxRadiusCumSum: seq[float] = newSeq[float](radii.len)
    diffFluxCDFs: seq[seq[float]] = newSeq[seq[float]](radii.len)
    diffRadiusSum = 0.0

  template toCdf(x: untyped): untyped =
    let baseline = x[0]
    let integral = x[^1]
    x.mapIt( (it - baseline) / (integral - baseline) )

  for iRad, radius in radii:
    # emRates is seq of radii of energies
    let emRate = emRates[iRad]
    var diffFlux = newSeq[float](emRate.len)
    var diffSum = 0.0
    var radiusCumSum = newSeq[float](energies.len)
    for iEnergy, energy in energies:
      diffFlux[iEnergy] = emRate[iEnergy] * (energy.float*energy.float) * radius*radius
      diffSum += diffFlux[iEnergy]
      radiusCumSum[iEnergy] = diffSum

    when false:
      # sanity checks for calc of differential flux & emission rate
      let df = toDf({ "diffFlux" : diffFlux, "energy" : energies.mapIt(it.float),
                      "emRate" : emRates[iRad] })
      ggplot(df, aes("energy", "diffFlux")) +
        geom_line() +
        scale_y_continuous() +
        ggsave("/tmp/diff_flux_vs_energy.pdf")
      ggplot(df, aes("energy", "emRate")) +
        geom_line() +
        scale_y_continuous() +
        ggsave("/tmp/emRate_vs_energy.pdf")
    diffRadiusSum += diffSum
    fluxRadiusCumSum[iRad] = diffRadiusSum
    diffFluxCDFs[iRad] = radiusCumSum.toCdf()
  let fluxRadiusCDF = fluxRadiusCumSum.toCdf()

  ## sample from random point and plot
  when false: # this is a sort of sanity check to check the sampling of 100_000 elements
    var es = newSeq[float]()
    var ems = newSeq[float]()
    var rs = newSeq[float]()
    var ts = newSeq[float]()
    var ps = newSeq[float]()

    for i in 0 ..< 100_000:
      let pos = getRandomPointFromSolarModel(centerSun, RadiusSun, emratesRadiusCumSum)
      let r = (pos - centerSun).length()
      let energyAx = getRandomEnergyFromSolarModel(
        pos, centerSun, RadiusSun, energies, emrates, diffFluxCDFs, "energy"
      )
      let em = getRandomEnergyFromSolarModel(
        pos, centerSun, RadiusSun, energies, emrates, diffFluxCDFs, "emissionRate"
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

  result = FullRaytraceSetup(
    centerVecs: centerVecs,
    expSetup: expSetup,
    detectorSetup: detectorSetup,
    energies: energies,
    fluxRadiusCDF: fluxRadiusCDF,
    diffFluxCDFs: diffFluxCDFs,
    flags: flags
  )

proc calculateFluxFractions(raytraceSetup: FullRaytraceSetup,
                            generatePlots = true, suffix = ""): seq[Axion] =
  ## In the following we will go over a number of points in the sun, whose location and
  ## energy will be biased by the emission rate and whose track will be through the CAST
  ## experimental setup from 2018 at VT3
  var axions = newSeq[Axion](NumberOfPointsSun)
  var axBuf = cast[ptr UncheckedArray[Axion]](axions[0].addr)
  echo "start"
  init(Weave)
  traceAxionWrapper(axBuf, NumberOfPointsSun,
                    raytraceSetup.centerVecs,
                    raytraceSetup.expSetup,
                    raytraceSetup.detectorSetup,
                    raytraceSetup.fluxRadiusCDF,
                    raytraceSetup.diffFluxCDFs,
                    raytraceSetup.energies,
                    raytraceSetup.flags)
  exit(Weave)

  if generatePlots:
    generateResultPlots(axions, raytraceSetup.detectorSetup.windowYear, suffix)
  result = axions

proc performAngularScan(angularScanMin, angularScanMax: float, numAngularScanPoints: int,
                        noPlots: bool,
                        flags: set[ConfigFlags]) =
  ## Performs a scan of the telescope efficiency (intended for the XMM Newton optics)
  ## under different angles.
  let (esKind, dkKind, skKind) = parseSetup()
  let angles = linspace(angularScanMin, angularScanMax, numAngularScanPoints)
  var fullSetup = initFullSetup(esKind,
                                dkKind,
                                skKind,
                                flags) # radiationCharacteristic = "axionRadiation::characteristic::sar"
  var fluxes = newSeq[float](numAngularScanPoints)
  for i, angle in angles:
    let suffix = &"angle_{angle:.2f}"
    # modify the angle of the telescope
    var expSetup = fullSetup.expSetup
    var tel = expSetup.telescope
    tel.telescope_turned_y = angle.Degree
    expSetup.telescope = tel
    fullSetup.expSetup = expSetup
    let axions = fullSetup.calculateFluxFractions(generatePlots = not noPlots, suffix = suffix)
    fluxes[i] = axions.filterIt(it.passed).mapIt(it.weights).sum()
  let maxFlux = fluxes.max
  fluxes.applyIt(it / maxFlux)
  var df = toDf({"Angle [deg]" : angles, "relative flux" : fluxes})
  echo df.pretty(-1)
  let dfXMM = readCsv("../resources/xmm_newton_angular_effective_area.csv", header = "#")
    .filter(f{`effectiveArea` > 0.0})
    .mutate(f{"Angle [deg]" ~ idx("angle[arcmin]") / 60.0},
            f{"effectiveArea" ~ `effectiveArea` / max(`effectiveArea`)})
  let dfMcXtrace = readCsv("../resources/McXtrace_angular_xmm.csv")
  df = bind_rows([("Nim", df), ("McXtrace", dfMcXtrace)], "Type")
  ggplot(df, aes("Angle [deg]", "relative flux", color = "Type")) +
    geom_point() +
    geom_line(data = dfXMM, aes = aes("Angle [deg]", "effectiveArea"), color = "#FF00FF") +
    ggtitle("Normalized total flux in scan of telescope angle. Solid line: XMM Newton 'theory'") +
    ggsave("../out/angular_scan_telescope_y.pdf", width = 800, height = 480)

proc main(
  ignoreDetWindow = false, ignoreGasAbs = false,
  ignoreConvProb = false, ignoreReflection = false, xrayTest = false,
  detectorInstall = false, magnet = false,
  angularScanMin = 0.0, angularScanMax = 0.0, numAngularScanPoints = 50,
  noPlots = false
         ) =
  # check if the `config.toml` file exists, otherwise recreate from the default
  if not fileExists("config.toml"):
    let cdata = readFile(sourceDir / "config_default.toml")
    writeFile("config.toml", cdata)

  var coldboreBlockedLength: float64
  coldboreBlockedLength = 0.0

  var flags: set[ConfigFlags]
  if ignoreDetWindow:  flags.incl cfIgnoreDetWindow
  if ignoreGasAbs:     flags.incl cfIgnoreGasAbs
  if ignoreConvProb:   flags.incl cfIgnoreConvProb
  if ignoreReflection: flags.incl cfIgnoreReflection
  if xrayTest: flags.incl cfXrayTest
  if magnet: flags.incl cfReadMagnetConfig
  if detectorInstall: flags.incl cfReadDetInstallConfig
  echo "Flags: ", flags

  if angularScanMin == angularScanMax:
    let (esKind, dkKind, skKind) = parseSetup()
    let fullSetup = initFullSetup(esKind,
                                  dkKind,
                                  skKind,
                                  flags) # radiationCharacteristic = "axionRadiation::characteristic::sar"
    discard fullSetup.calculateFluxFractions(generatePlots = not noPlots)
  else:
    # perform a scan of the angular rotation of the telescope
    performAngularScan(angularScanMin, angularScanMax, numAngularScanPoints, noPlots, flags)

when isMainModule:
  dispatch main
