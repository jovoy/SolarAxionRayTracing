import glm/vec
import math
import random
# kompilieren und ausführen: nim c -r ael.nim

import nim-plotly-master/src/plotly
import nim-plotly-master/src/plotly/chroma
import nim-plotly-master/src/plotly/names
import random
import sequtils

#nim c -r fig9_heatmap.nim

let
  # The GL heatmap is also supported as HeatMapGL
  d = Trace[float32](mode: PlotMode.Lines, `type`: PlotType.HeatMap)

d.colormap = ColorMap.Viridis
# fill data for colormap with random values. The data needs to be supplied
# as a nested seq.
d.zs = newSeqWith(28, newSeq[float32](28))
for x in 0 ..< 28:
  for y in 0 ..< 28:
    d.zs[x][y] = random(1.0)
let
  layout = Layout(title: "Heatmap example", width: 800, height: 800,
                  xaxis: Axis(title: "A heatmap x-axis"),
                  yaxis: Axis(title: "y-axis too"), autosize: false)
  p = Plot[float32](layout: layout, traces: @[d])
echo p.save()
#p.show()




##################rayTracer###############################
# ToDo: Picture heatmap
# put in XRT effeciency



################################
# VARIABLES from rayTracer.h
var RAYTRACER_DISTANCE_SUN_EARTH =  1.5e14  #mm
var RAYTRACER_RADIUS_SUN = 6.9e11 #mm
var RAYTRACER_RADIUS_COLDBORE = 21.5 #mm
var RAYTRACER_LENGTH_COLDBORE = 9756.0 #mm half B field to end of CB
var RAYTRACER_LENGTH_COLDBORE_9T = 9260.0 #mm half B field to half B field
var RAYTRACER_LENGTH_PIPE_CB_VT3 = 2571.5 #mm should stay the same #from beam pipe drawings
var RAYTRACER_RADIUS_PIPE_CB_VT3 = 33.6 #mm from drawing #30.0 #mm smallest aperture between end of CB and VT4
var RAYTRACER_LENGTH_PIPE_VT3_XRT = 264.7 #mm from picture #198.2 #mm from XRT drawing
var RAYTRACER_RADIUS_PIPE_VT3_XRT = 25.0 #mm from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm
var RAYTRACER_FOCAL_LENGTH_XRT = 1500.0 #mm is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 #mm was the Telescope of 2014 (MPE XRT) also: Aperatur changed
var RAYTRACER_DISTANCE_AXIS_CB_AXIS_XRT = 58.44 #mm from XRT drawing #no change, because don't know
var RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW = -10.0 #mm #no change, because don't know
var numberOfPointsEndOfCB = 1000
var numberOfPointsSun = 1000

## Chipregions#####

var CHIPREGIONS_CHIP_X_MIN = 0.0
var CHIPREGIONS_CHIP_X_MAX = 14.0
var CHIPREGIONS_CHIP_Y_MIN = 0.0
var CHIPREGIONS_CHIP_Y_MAX = 14.0

var CHIPREGIONS_CHIP_CENTER_X = 7.0
var CHIPREGIONS_CHIP_CENTER_Y = 7.0

var CHIPREGIONS_GOLD_X_MIN = 4.5
var CHIPREGIONS_GOLD_X_MAX = 9.5
var CHIPREGIONS_GOLD_Y_MIN = 4.5
var CHIPREGIONS_GOLD_Y_MAX = 9.5

var CHIPREGIONS_SILVER_RADIUS_MAX = 4.5
var CHIPREGIONS_BRONZE_RADIUS_MAX = 5.5

###### from rootConfig.h
var ROOTCONFIG_CANVAS_WIDTH = 1600.0
var ROOTCONFIG_CANVAS_HEIGHT = 1000.0
var ROOTCONFIG_SQUARECANVAS_WIDTH = 720.0
var ROOTCONFIG_SQUARECANVAS_HEIGHT = 600.0
var ROOTCONFIG_CANVAS_MARGIN_LEFT = 0.15
var ROOTCONFIG_CANVAS_MARGIN_RIGHT = 0.05
var ROOTCONFIG_CANVAS_MARGIN_BOTTOM = 0.15
var ROOTCONFIG_CANVAS_MARGIN_TOP = 0.1
var ROOTCONFIG_SQUARECANVAS_MARGIN_LEFT = 0.15
var ROOTCONFIG_SQUARECANVAS_MARGIN_RIGHT = 0.225
var ROOTCONFIG_SQUARECANVAS_MARGIN_BOTTOM = 0.15
var ROOTCONFIG_SQUARECANVAS_MARGIN_TOP = 0.1
var ROOTCONFIG_FONT_TIMESNEWROMAN = 132
var ROOTCONFIG_FONT_ARIAL = 42
################################

proc rayTracer(numberOfPointsEndOfCB : int, numberOfPointsSun : int) : int= 1
    
   
    
var fluxFractionGold = 0.0 #dies muss eine globale var sein
proc getFluxFractionGold(): float64 =
  #fluxFractionGold = 0.0
  result = fluxFractionGold

var fluxFractionSilver = 0.0 
proc getFluxFractionSilver(): float64 =
  result = fluxFractionSilver

var fluxFractionBronze = 0.0 
proc getFluxFractionBronze(): float64 =
  result = fluxFractionBronze
                        
                  # der Typ ist noch nicht definiert <-- machen wir später, wenn wir Objekte bauen.Mit String ist auch ok, nur nicht so fancy
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

#echo getFluxFraction("region: silver")
#echo getFluxFraction("region: gold")
#echo getFluxFraction("region: goldPlusSilver")
                            # needs to be object
proc calculateFluxFractions(axionRadiationCharacteristic: string, 
  xrtTransmissionAt10Arcmin :float64,
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
  # ok, wir definieren das im global scope hier in der Datei und müssen das dann irgendwann umziehen in die Hauptdatei

  var misalignmentDetector = vec3(0.0)
  misalignmentDetector[0] = detectorMisalignmentX
  misalignmentDetector[1] = detectorMisalignmentY
  misalignmentDetector[2] = 0

  var centerSun = vec3(0.0)
  centerSun[0] = 0
  centerSun[1] = 0
  centerSun[2] = RAYTRACER_DISTANCE_SUN_EARTH
  
  centerSun = centerSun + misalignmentSun
  #echo centerSun

  var radiusSun = RAYTRACER_RADIUS_SUN

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

  var radiusCB = RAYTRACER_RADIUS_COLDBORE

  var centerExitPipeCBVT3 = vec3(0.0)
  centerExitPipeCBVT3[0] = 0
  centerExitPipeCBVT3[1] = 0
  centerExitPipeCBVT3[2] = RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3

  var radiusPipeCBVT3 = RAYTRACER_RADIUS_PIPE_CB_VT3

  var centerExitPipeVT3XRT = vec3(0.0)
  centerExitPipeVT3XRT[0] = 0
  centerExitPipeVT3XRT[1] = 0
  centerExitPipeVT3XRT[2] = RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT

  var radiusPipeVT3XRT = RAYTRACER_RADIUS_PIPE_VT3_XRT

  var distanceCBAxisXRTAxis = RAYTRACER_DISTANCE_AXIS_CB_AXIS_XRT

  var integralNormalisation = 0.0
  var integralTotal = 0.0
  var integralDetector = 0.0
  var integralBronze = 0.0
  var integralSilver = 0.0
  var integralGold = 0.0

  var misalignment = (sunMisalignmentH != 0.0) or (sunMisalignmentV != 0.0) or (detectorMisalignmentX != 0.0) or (detectorMisalignmentY != 0.0)
  # am Anfang fehlt noch eteas
  var plotName = ""
  plotname = plotname & axionRadiationCharacteristic #hier müssen wir wieder was überlegen, da das Objekt in den Klammern nochnichtnachgebaut ist
  plotname = plotname & "-txrtXX-" 
  var xrtTransmission_path_int = xrtTransmissionAt10Arcmin * 100.0 + 0.5 #//das ist eine übergebene Variable -> ok
  echo xrtTransmission_path_int
  plotname = plotname & $xrtTransmission_path_int

  if misalignment == true:
    plotName = plotname & "-" & $(sunMisalignmentH) & "'"
    plotName = plotname & "-" & $(sunMisalignmentV) & "'"
    plotName = plotname & "-" & $(detectorMisalignmentX) & "mm"
    plotName = plotname & "-" & $(detectorMisalignmentY) & "mm"

  plotName = plotname & "-axion-image"

  #canvas graphics area dimensions but for root part missing
  var w = ROOTCONFIG_SQUARECANVAS_WIDTH
  var h = ROOTCONFIG_SQUARECANVAS_HEIGHT
  var l = ROOTCONFIG_SQUARECANVAS_MARGIN_LEFT
  var r = ROOTCONFIG_SQUARECANVAS_MARGIN_RIGHT
  var b = ROOTCONFIG_SQUARECANVAS_MARGIN_BOTTOM
  var t = ROOTCONFIG_SQUARECANVAS_MARGIN_TOP
  ##

  ######### some functions we're gonna need later##########
  # Now let's define a function to get us a random point on a disk, where we later can put the center of the
  # detector as center and the radius of the detector as radius
  
  proc getRandomPointOnDisk(center: Vec3, radius:float64) : Vec3 =
    var x = 0.0
    var y = 0.0
    var r = radius * sqrt(random(1.0)) #Zufallszahl zwischen 0 und 1 (hoffe), die Wurzel nach C++-Vorlage(wtf)
    # _randomGEnerator -> Circle(x,y,r)  #wird nicht weiter verwendet (?)
    var angle = 360 * random(1.0) #random angle
    x = cos(angle) * r
    y = sin(angle) * r
    var vector = vec3(x,y,0.0)
    vector = vector + center
    result = vector 

  proc getRandomPointFromSolarModel(center : Vec3, radius : float64 ) : Vec3 =
    var x = 0.0
    var y = 0.0
    var z = 0.0
    var r = radius * 1e-1 * random(1.0)
    var angle1 = 360 * random(1.0)
    var angle2 = 180 * random(1.0)
    x = cos(angle1) * sin(angle2) * r
    y = sin(angle1) * sin(angle2) * r
    z = cos(angle2) * r
    var vector = vec3(x,y,z)
    vector = vector + center
    result = vector

  echo getRandomPointFromSolarModel(centerSun,radiusSun) #Bedenke: der Ursprung dieser Koordinaten ist die Erde
  

  proc lineIntersectsCircle(point_1 : Vec3, point_2 : Vec3, center : Vec3, radius : float64, intersect : Vec3) : bool =
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
    var intersect_2 = vec3(0.0)
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

  proc telescopeTransmission(angle : float64,transmissionat10Arcmin :float64) : float64 =   # this has to be changed since it is the old XRT
    var b : float64
    b = 1.0
    var a : float64
    a = (1.0 - transmissionAt10Arcmin) / (0.0 - degToRad(10.0/60.0) )
    var t : float64
    t = a*angle + b
    if (t<0.0):
      t = 0.0
      return t 

  ############done with the functions, let's use them############

  for iExitCB in 0..<numberOfPointsEndOfCB:
    #echo iExitCB
    var pointExitCBMagneticField = vec3(0.0)
    pointExitCBMagneticField = getRandomPointOnDisk(centerExitCBMagneticField, radiusCB)
    for iSun in 0..<numberOfPointsSun:
        # echo numberOfPointsSun
        integralNormalisation = integralNormalisation + 1
        var pointInSun = vec3(0.0)
        case  axionRadiationCharacteristic
        of "axionRadiation::characteristic::sar":
          pointInSun = getRandomPointFromSolarModel(centerSun,radiusSun)
        of "axionRadiation::characteristic::def":
          echo "Error: Default radiation characteristic not implemented"
        else:
          echo "Error: Unknown axion radiation characteristic"

        var intersect = vec3(0.0)
        var pathCB : float64
        var intersectsEntranceCB : bool
        intersectsEntranceCB = lineIntersectsCircle(pointInSun, pointExitCBMagneticField, centerEntranceCB, radiusCB, intersect)
        var intersectsCB : bool
        intersectsCB = false

        if (not intersectsEntranceCB):
            intersectsCB = lineIntersectsCylinderOnce(pointInSun,pointExitCBMagneticField,centerEntranceCB,centerExitCBMagneticField,radiusCB,intersect)

        if (not intersectsEntranceCB and not intersectsCB): continue

        pathCB = pointExitCBMagneticField[2] - intersect[2]
        var pointExitCB = vec3(0.0)
        
        if (not lineIntersectsCircle(pointInSun,pointExitCBMagneticField,centerExitCB,radiusCB,pointExitCB)): continue

        var pointExitPipeCBVT3 =vec3(0.0)

        if (not lineIntersectsCircle(pointExitCBMagneticField,pointExitCB,centerExitPipeCBVT3,radiusPipeCBVT3,pointExitPipeCBVT3)) : continue

        var pointExitPipeVT3XRT = vec3(0.0)

        if (not lineIntersectsCircle(pointExitCB,pointExitPipeCBVT3,centerExitPipeVT3XRT,radiusPipeVT3XRT,pointExitPipeVT3XRT)) : continue

        var vectorBeforeXRT = vec3(0.0)
        vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB

        ########von CB zum XRT

        var pointEntranceXRT = vec3(0.0)
        pointEntranceXRT[0] = pointExitPipeVT3XRT[0]
        pointEntranceXRT[1] = pointExitPipeVT3XRT[1] + distanceCBAxisXRTAxis
        pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

        var angle : float64
        angle = vectorBeforeXRT[1]  # Here we want to adress theta, the polar angle, which should be the second entrance of the vector
        var r_x : float64
        r_x = pointEntranceXRT[0]
        var r_y : float64
        r_y = pointEntranceXRT[1]
        var theta_x : float64
        theta_x = vectorBeforeXRT[0] / vectorBeforeXRT[2]
        var theta_y : float64
        theta_y = vectorBeforeXRT[1] / vectorBeforeXRT[2]
        var theta_x_prime : float64
        theta_x_prime= theta_x - ( r_x / RAYTRACER_FOCAL_LENGTH_XRT) 
        var theta_y_prime : float64
        theta_y_prime= theta_y - ( r_y / RAYTRACER_FOCAL_LENGTH_XRT)

        var vectorAfterXRT =vec3(0.0)
        vectorAfterXRT[0] = sin(theta_x_prime) * 100.0
        vectorAfterXRT[1] = sin(theta_y_prime) * 100.0
        vectorAfterXRT[2] = 100.0
        echo "now" 
        echo vectorAfterXRT

        var centerDetectorWindow = vec3(0.0)
        centerDetectorWindow[0] = 0.0
        centerDetectorWindow[1] = 0.0
        centerDetectorWindow[2] = pointEntranceXRT[2] + RAYTRACER_FOCAL_LENGTH_XRT + RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW

        var lambda_0 : float64
        lambda_0 = ( centerDetectorWindow[2] - pointEntranceXRT[2] ) / vectorAfterXRT[2]
        var pointDetectorWindow = vec3(0.0)
        pointDetectorWindow = pointEntranceXRT + lambda_0 * vectorAfterXRT
        pointDetectorWindow = pointDetectorWindow - misalignmentDetector

        var weight : float64
        weight = ( telescopeTransmission(angle,xrtTransmissionAt10Arcmin) * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) )
        
        integralTotal = integralTotal + weight

        ###detector COS has (0/0) at the bottom left corner of the chip

        pointDetectorWindow[0] = pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X
        pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y

        var gold : bool
        gold = ( (pointDetectorWindow[0] >= CHIPREGIONS_GOLD_X_MIN) and (pointDetectorWindow[0] <= CHIPREGIONS_GOLD_X_MAX) and (pointDetectorWindow[1] >= CHIPREGIONS_GOLD_Y_MIN) and (pointDetectorWindow[1] <= CHIPREGIONS_GOLD_Y_MAX) )

        var r_xy : float64
        r_xy = sqrt( ( (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X) * (pointDetectorWindow[0] - CHIPREGIONS_CHIP_CENTER_X) ) + ( (pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y) * (pointDetectorWindow[1] - CHIPREGIONS_CHIP_CENTER_Y) ) )

        var silver : bool
        silver = (r_xy <= CHIPREGIONS_SILVER_RADIUS_MAX ) and not gold

        var bronze : bool
        bronze = not gold and not silver and ( r_xy <= CHIPREGIONS_BRONZE_RADIUS_MAX )

        var withinWindow : bool
        withinWindow = r_xy < detectorWindowAperture/2

        var detector : bool
        detector = ( (pointDetectorWindow[0] >= CHIPREGIONS_CHIP_X_MIN) and (pointDetectorWindow[0] <= CHIPREGIONS_CHIP_X_MAX) and (pointDetectorWindow[1] >= CHIPREGIONS_CHIP_Y_MIN) and (pointDetectorWindow[1] <= CHIPREGIONS_CHIP_Y_MAX) )

        if(gold and withinWindow): integralGold = integralGold + weight
        if(silver and withinWindow): integralSilver = integralSilver + weight
        if(bronze and withinWindow): integralBronze = integralBronze + weight
        if(detector and withinWindow): integralDetector = integralDetector + weight
        
        #if(withinWindow){  image->Fill(pointDetectorWindow[0],pointDetectorWindow[1],weight)









       


########################### aEL main ############################

#configuration
var showExpectedLimitOnly : bool 
showExpectedLimitOnly = true
var plotExpectedLimitOnly : bool 
plotExpectedLimitOnly = false

#axion specific stuff
var radiationCharacteristic : string ##axionRadiation::characteristic radiationCharacteristic(axionRadiation::characteristic::sar);
radiationCharacteristic = "axionRadiation::characteristic::sar"
var coldboreBlockedLength : float64 
coldboreBlockedLength = 0.0

#Xrt specific stuff
var perfectXRT : bool 
perfectXRT = false #if true transmission of 100 % is assumed for all energies else interpolated data from transmission measurements is used
var xrtTransmissionAt10Arcmin : float64 
xrtTransmissionAt10Arcmin = 0.7 #relative transmission for x-rays at 10' angle compared to parallel beams #need to be changed?

#detector specific stuff
#detector::detectorVersion detectorVersion(detector::detectorVersion::mk2);
#detector::windowMaterial detectorWindowMaterial(detector::windowMaterial::Mylar);
var detectorWindowThickness : float64 
detectorWindowThickness = 2000.0 #nm
var detectorWindowAreaTransparency : float64 
detectorWindowAreaTransparency = 0.826 #Mylar 3x3: 82.6 %, SiN Quad-Ribs 500 µm: 73.5 %
var detectorWindowAluminiumLayerThickness : float64 
detectorWindowAluminiumLayerThickness = 40.0 #nm Mylar: 40 nm, SiN: 20 nm
#detector::windowMaterial differentialWindowMaterial(detector::windowMaterial::Mylar);
var differentialWindowThickness :float64
differentialWindowThickness = 900.0 #nm
var detectorWindowAperture : float64 
detectorWindowAperture = 14.0 #mm
var smearEnergy : bool 
smearEnergy = false  #true
        
#analysis specific stuff
var useSystematicErrors : bool 
useSystematicErrors = false   #true
var useStatisticalErrors : bool 
useStatisticalErrors = false
var createBackgroundRatePlots : bool 
createBackgroundRatePlots = true


echo calculateFluxFractions(radiationCharacteristic, xrtTransmissionAt10Arcmin, detectorWindowAperture, 0.0, 0.0, 0.0, 0.0, coldboreBlockedLength) # radiationCharacteristic = "axionRadiation::characteristic::sar"

#echo calculateFluxFractions(,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
                          #hier muss der richtige String rein

#type
#    vector3 = ref object of Vec3

## things changed##
# XRT Focal length 1600.0 -> 1500.0
