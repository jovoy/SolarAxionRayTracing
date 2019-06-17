import glm/vec
import math
import random
# kompilieren und ausführen: nim cpp -r aEL.nim

import nim-plotly-master/src/plotly
import nim-plotly-master/src/plotly/chroma
import nim-plotly-master/src/plotly/names
import random
import sequtils
import nimhdf5-master/src/nimhdf5
import nimhdf5-master/src/nimhdf5/H5nimtypes
import nimhdf5-master/src/nimhdf5/hdf5_wrapper

#nim cpp -r fig9_heatmap.nim




##################rayTracer###############################
# ToDo: Picture heatmap

# put in XRT effeciency



################################
# VARIABLES from rayTracer.h
const RAYTRACER_DISTANCE_SUN_EARTH =  1.5e14  #mm
const RAYTRACER_RADIUS_SUN = 6.9e11 #mm
const RAYTRACER_RADIUS_COLDBORE = 21.5 #mm
const RAYTRACER_LENGTH_COLDBORE = 9756.0 #mm half B field to end of CB
const RAYTRACER_LENGTH_COLDBORE_9T = 9260.0 #mm half B field to half B field
const RAYTRACER_LENGTH_PIPE_CB_VT3 = 2571.5 #mm should stay the same #from beam pipe drawings
const RAYTRACER_RADIUS_PIPE_CB_VT3 = 30.0#33.6 #mm from drawing #30.0 #mm smallest aperture between end of CB and VT4
const RAYTRACER_LENGTH_PIPE_VT3_XRT = 198.2#264.7 #mm from picture #198.2 #mm from XRT drawing
const RAYTRACER_RADIUS_PIPE_VT3_XRT = 35.0#25.0 #mm from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm
const RAYTRACER_FOCAL_LENGTH_XRT = 1600.0#1500.0 #mm is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 #mm was the Telescope of 2014 (MPE XRT) also: Aperatur changed
const RAYTRACER_DISTANCE_AXIS_CB_AXIS_XRT = 58.44 #mm from XRT drawing #no change, because don't know
const RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW = -10.0 #mm #no change, because don't know
const numberOfPointsEndOfCB = 1000
const numberOfPointsSun = 1000
var xrtTransmissionAt10Arcmin : float64 
xrtTransmissionAt10Arcmin = 0.7 #relative transmission for x-rays at 10' angle compared to parallel beams #need to be changed?
## Chipregions#####

const CHIPREGIONS_CHIP_X_MIN = 0.0
const CHIPREGIONS_CHIP_X_MAX = 14.0
const CHIPREGIONS_CHIP_Y_MIN = 0.0
const CHIPREGIONS_CHIP_Y_MAX = 14.0

const CHIPREGIONS_CHIP_CENTER_X = 7.0
const CHIPREGIONS_CHIP_CENTER_Y = 7.0

const CHIPREGIONS_GOLD_X_MIN = 4.5
const CHIPREGIONS_GOLD_X_MAX = 9.5
const CHIPREGIONS_GOLD_Y_MIN = 4.5
const CHIPREGIONS_GOLD_Y_MAX = 9.5

const CHIPREGIONS_SILVER_RADIUS_MAX = 4.5
const CHIPREGIONS_BRONZE_RADIUS_MAX = 5.5

###### from rootConfig.h
const ROOTCONFIG_CANVAS_WIDTH = 1600.0
const ROOTCONFIG_CANVAS_HEIGHT = 1000.0
const ROOTCONFIG_SQUARECANVAS_WIDTH = 720.0
const ROOTCONFIG_SQUARECANVAS_HEIGHT = 600.0
const ROOTCONFIG_CANVAS_MARGIN_LEFT = 0.15
const ROOTCONFIG_CANVAS_MARGIN_RIGHT = 0.05
const ROOTCONFIG_CANVAS_MARGIN_BOTTOM = 0.15
const ROOTCONFIG_CANVAS_MARGIN_TOP = 0.1
const ROOTCONFIG_SQUARECANVAS_MARGIN_LEFT = 0.15
const ROOTCONFIG_SQUARECANVAS_MARGIN_RIGHT = 0.225
const ROOTCONFIG_SQUARECANVAS_MARGIN_BOTTOM = 0.15
const ROOTCONFIG_SQUARECANVAS_MARGIN_TOP = 0.1
const ROOTCONFIG_FONT_TIMESNEWROMAN = 132
const ROOTCONFIG_FONT_ARIAL = 42
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

var fluxFractionDetector = 0.0 
proc getFluxFractionDetector(): float64 =
  result = fluxFractionDetector

var fluxFractionTotal = 0.0 
proc getFluxFractionTotal(): float64 =
  result = fluxFractionTotal
                        
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
  # , vecdata_x : array[0..1000, float], vecdata_y : array[0..1000,float]

  proc prepareheatmap( numberofrows : int, numberofcolumns : int, start_x : float, stop_x : float, start_y : float, stop_y : float, data_X :seq, data_Y : seq, weights : seq) : any =
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
    var heatmaptable = newSeqWith(numberofrows, newSeq[float](numberofcolumns))
    for i, value in data_X:
      var coord_X = floor((data_X[i] - start_x) / stepsize_X)
      var coord_Y = floor((data_Y[i] - start_y) / stepsize_Y)
      if coord_X > 0:
        if coord_Y > 0:
          if coord_X < float(numberofrows):
            if coord_Y < float(numberofcolumns):
              heatmaptable[int(coord_X)][int(coord_Y)] = heatmaptable[int(coord_X)][int(coord_Y)] + 1
    result = heatmaptable


  proc drawfancydiagrams(diagramtitle : string, objectstodraw : any, width : int) : float =
    let
      # The GL heatmap is also supported as HeatMapGL
      d = Trace[float32](mode: PlotMode.Lines, `type`: PlotType.HeatMap)
    
    d.colormap = ColorMap.Viridis
    # fill data for colormap with random values. The data needs to be supplied
    # as a nested seq.
    #
    #
    #
    #
    #
    d.zs = newSeqWith(width, newSeq[float32](width))
    for x in 0 ..< width:
      for y in 0 ..< width:
        if x < width:
          if y < width:
            if x > 0:
              if y > 0:
                d.zs[x][y] = objectstodraw[x][y] 
  
    const 
      y = @[float32(CHIPREGIONS_GOLD_Y_MIN * 40.0 / 14.0), float32(CHIPREGIONS_GOLD_Y_MIN * 40.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 40.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 40.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MIN * 40.0 / 14.0)]
      x = @[float32(CHIPREGIONS_GOLD_X_MIN * 40.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 40.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 40.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 40.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 40.0 / 14.0)]
    let
      d4 = Trace[float32](mode: PlotMode.LinesMarkers, `type`: PlotType.ScatterGL, ys: y, xs : x) 
      
                     
    let
      layout = Layout(title: diagramtitle, width: 800, height: 800,
                      xaxis: Axis(title: "x-axis 40 = 14 mm"),
                      yaxis: Axis(title: "y-axis 40 = 14 mm"), autosize: false)
      p = Plot[float32](layout: layout, traces: @[d, d4])
    echo p.save()
    p.show()
    var b = 2.2
    result = b


  #var ara1 = @[31.0,23.0]
  #var ara2 = @[21.0,12.0]
  #var dataX = @[1.1,2.1,1.3,2.3,2.4,1.3,7.6,8.9,8.9,8.9,8.9]
  #var dataY = @[4.1,3.1,1.3,2.3,2.4,1.3,4.5,1.8,1.6,1.6,1.6]
  #var heatmaptable = prepareheatmap(10,10,10.0,0.0,10.0,0.0,dataX,dataY)
  #echo heatmaptable
  #echo drawfancydiagrams("Diagramtitel", heatmaptable, 10)

  var array1 = [[1,2],[3,4]]
  echo array1[1][1]
  echo array1[1][0]
  echo array1[0][1]
   
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
  

  proc lineIntersectsCircle(point_1 : Vec3, point_2 : Vec3, center : Vec3, radius : float64, intersect : Vec3) : bool = # probably still some error (lambda1 -> infinity)
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
  
  var pointdataX = @[0.0] 
  var pointdataY = @[0.0]
  var weights = @[0.0]
  
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
        #echo pointInSun
        var intersect = vec3(0.0)
        var pathCB : float64
        var intersectsEntranceCB : bool
        intersectsEntranceCB = lineIntersectsCircle(pointInSun, pointExitCBMagneticField, centerEntranceCB, radiusCB, intersect)
        var intersectsCB : bool
        intersectsCB = false
        
        if (not intersectsEntranceCB):
            intersectsCB = lineIntersectsCylinderOnce(pointInSun,pointExitCBMagneticField,centerEntranceCB,centerExitCBMagneticField,radiusCB,intersect)
          
        if (not intersectsEntranceCB and not intersectsCB): continue

        intersect = pointInSun + ((centerEntranceCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
        
        pathCB = pointExitCBMagneticField[2] - intersect[2]
        var pointExitCB = vec3(0.0)
        
        if (not lineIntersectsCircle(pointInSun,pointExitCBMagneticField,centerExitCB,radiusCB,pointExitCB)): continue
        
        pointExitCB = pointInSun + ((centerExitCB[2] - pointInSun[2]) / (pointExitCBMagneticField - pointInSun)[2]) * (pointExitCBMagneticField - pointInSun)
        var pointExitPipeCBVT3 =vec3(0.0)
        
        if (not lineIntersectsCircle(pointExitCBMagneticField,pointExitCB,centerExitPipeCBVT3,radiusPipeCBVT3,pointExitPipeCBVT3)) : continue
        
        pointExitPipeCBVT3 = pointExitCBMagneticField + ((centerExitPipeCBVT3[2] - pointExitCBMagneticField[2]) / (pointExitCB - pointExitCBMagneticField)[2]) * (pointExitCB - pointExitCBMagneticField)
        var pointExitPipeVT3XRT = vec3(0.0)
        
        if (not lineIntersectsCircle(pointExitCB,pointExitPipeCBVT3,centerExitPipeVT3XRT,radiusPipeVT3XRT,pointExitPipeVT3XRT)) : continue

        pointExitPipeVT3XRT = pointExitCB + ((centerExitPipeVT3XRT[2] - pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) * (pointExitPipeCBVT3 - pointExitCB)
        
        #pointExitPipeVT3XRT = pointInSun - pointExitCBMagneticField 
        #pointExitPipeVT3XRT = pointExitPipeVT3XRT/ sqrt(pointExitPipeVT3XRT[0]*pointExitPipeVT3XRT[0] + pointExitPipeVT3XRT[1]*pointExitPipeVT3XRT[1] + pointExitPipeVT3XRT[2]*pointExitPipeVT3XRT[2])
        #pointExitPipeVT3XRT = pointExitCBMagneticField + pointExitPipeVT3XRT * (RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT)

        var vectorBeforeXRT = vec3(0.0)
        vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB
        #vectorBeforeXRT[0] = 0.4
        #vectorBeforeXRT[1] = -0.2
        #vectorBeforeXRT[2] = 2769.7
        
        ###################von CB zum XRT#######################

        var pointEntranceXRT = vec3(0.0)
        pointEntranceXRT[0] = pointExitPipeVT3XRT[0]
        pointEntranceXRT[1] = pointExitPipeVT3XRT[1] + distanceCBAxisXRTAxis
        pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

        var angle : float64
        angle = (arccos(vectorBeforeXRT[2]/sqrt(vectorBeforeXRT[0]*vectorBeforeXRT[0]+vectorBeforeXRT[1]*vectorBeforeXRT[1]+vectorBeforeXRT[2]*vectorBeforeXRT[2])))  # Here we want to adress theta, the polar angle, which should be the second entrance of the vector
        var r_x : float64#vectorBeforeXRT[1] #
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
        
        var vectorAfterXRT = vec3(0.0)
        vectorAfterXRT[0] = sin(theta_x_prime) * 100.0
        vectorAfterXRT[1] = sin(theta_y_prime) * 100.0
        vectorAfterXRT[2] = 100.0
        
        var vectorAfterXRTPolar = vec3(0.0)  #(r,theta,phi)
        vectorAfterXRTPolar[0] = sqrt(vectorAfterXRT[0]*vectorAfterXRT[0]+vectorAfterXRT[1]*vectorAfterXRT[1]+vectorAfterXRT[2]*vectorAfterXRT[2])
        vectorAfterXRTPolar[1] = radToDeg(arccos(vectorAfterXRT[0]/vectorAfterXRTPolar[0]))
        vectorAfterXRTPolar[2] = radToDeg(arctan2(vectorAfterXRT[2],vectorAfterXRT[1]))

        vectorAfterXRTPolar[1] = 90.0 - vectorAfterXRTPolar[1] #this is the pitch angle
        vectorAfterXRTPolar[2] = 90.0 - vectorAfterXRTPolar[2] #this is the yaw angle

        #echo "now" 
        #echo angle
        #echo vectorAfterXRTPolar


        var centerDetectorWindow = vec3(0.0)
        centerDetectorWindow[0] = 0.0
        centerDetectorWindow[1] = 0.0
        centerDetectorWindow[2] = pointEntranceXRT[2] + RAYTRACER_FOCAL_LENGTH_XRT + RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW

        var lambda_0 : float64
        lambda_0 = ( centerDetectorWindow[2] - pointEntranceXRT[2] ) / vectorAfterXRT[2]
        var pointDetectorWindow = vec3(0.0)
        pointDetectorWindow = pointEntranceXRT + lambda_0 * vectorAfterXRT
        pointDetectorWindow = pointDetectorWindow - misalignmentDetector

        #echo telescopeTransmission(angle,xrtTransmissionAt10Arcmin)
        #echo (1.0 - xrtTransmissionAt10Arcmin) / (0.0 - degToRad(10.0/60.0) )
        var weight : float64
        weight = ( telescopeTransmission(angle,xrtTransmissionAt10Arcmin) * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) )
        #echo weight
        integralTotal = integralTotal + weight
        
        ###detector COS has (0/0) at the bottom left corner of the chip

        pointDetectorWindow[0] = pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X
        pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y
        
        pointdataX.add(pointDetectorWindow[0])
        pointdataY.add(pointDetectorWindow[1])
        weights.add(weight)


        #echo pointDetectorWindow

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
  #echo pointdataX
  #echo pointdataY    
  var heatmaptable2 = prepareheatmap(40,40,5.0,9.0,5.4,9.4,pointdataX,pointdataY,weights)
  #echo heatmaptable2
  echo drawfancydiagrams("Diagramtitel", heatmaptable2, 40)

  fluxFractionTotal    = integralTotal    / integralNormalisation
  fluxFractionDetector = integralDetector / integralNormalisation
  fluxFractionBronze   = integralBronze   / integralNormalisation
  fluxFractionSilver   = integralSilver   / integralNormalisation
  fluxFractionGold     = integralGold     / integralNormalisation

  echo "Flux fraction for the gold region:"
  echo getFluxFraction("region: gold")

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
#var xrtTransmissionAt10Arcmin : float64 
#xrtTransmissionAt10Arcmin = 0.7 #relative transmission for x-rays at 10' angle compared to parallel beams #need to be changed?

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

#########################efficiencyCalculator#####################

#constants#

const EFFICIENCYCALCULATOR_NUMBER_RANGES = 8
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE0 = "calibration-cdl-apr2014-C-EPIC-0.6kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE0_LOW = 0.15
const EFFICIENCYCALCULATOR_ENERGY_RANGE0_HIGH = 0.4
const EFFICIENCYCALCULATOR_CHARGE_RANGE0_MIN = 0.0
const EFFICIENCYCALCULATOR_CHARGE_RANGE0_MAX = 5.0e4
const EFFICIENCYCALCULATOR_RMSY_RANGE0_MIN = 0.1
const EFFICIENCYCALCULATOR_RMSY_RANGE0_MAX = 20.0
const EFFICIENCYCALCULATOR_LENGTH_RANGE0_MAX = 6.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE1 = "calibration-cdl-apr2014-Cu-EPIC-0.9kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE1_LOW = 0.4
const EFFICIENCYCALCULATOR_ENERGY_RANGE1_HIGH = 0.7
const EFFICIENCYCALCULATOR_CHARGE_RANGE1_MIN = 3.0e4
const EFFICIENCYCALCULATOR_CHARGE_RANGE1_MAX = 8.0e4
const EFFICIENCYCALCULATOR_RMSY_RANGE1_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE1_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE1_MAX = 6.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE2 = "calibration-cdl-apr2014-Cu-EPIC-2kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE2_LOW = 0.7
const EFFICIENCYCALCULATOR_ENERGY_RANGE2_HIGH = 1.2
const EFFICIENCYCALCULATOR_CHARGE_RANGE2_MIN = 7.0e4
const EFFICIENCYCALCULATOR_CHARGE_RANGE2_MAX = 1.3e5
const EFFICIENCYCALCULATOR_RMSY_RANGE2_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE2_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE2_MAX = 7.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE3 = "calibration-cdl-apr2014-Al-Al-4kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE3_LOW = 1.2
const EFFICIENCYCALCULATOR_ENERGY_RANGE3_HIGH = 2.1
const EFFICIENCYCALCULATOR_CHARGE_RANGE3_MIN = 9.0e4
const EFFICIENCYCALCULATOR_CHARGE_RANGE3_MAX = 2.1e5
const EFFICIENCYCALCULATOR_RMSY_RANGE3_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE3_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE3_MAX = 7.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE4 = "calibration-cdl-apr2014-Ag-Ag-6kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE4_LOW = 2.1
const EFFICIENCYCALCULATOR_ENERGY_RANGE4_HIGH = 3.2
const EFFICIENCYCALCULATOR_CHARGE_RANGE4_MIN = 2.0e5
const EFFICIENCYCALCULATOR_CHARGE_RANGE4_MAX = 4.0e5
const EFFICIENCYCALCULATOR_RMSY_RANGE4_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE4_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE4_MAX = 7.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE5 = "calibration-cdl-apr2014-Ti-Ti-9kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE5_LOW = 3.2
const EFFICIENCYCALCULATOR_ENERGY_RANGE5_HIGH = 4.9
const EFFICIENCYCALCULATOR_CHARGE_RANGE5_MIN = 2.9e5
const EFFICIENCYCALCULATOR_CHARGE_RANGE5_MAX = 5.5e5
const EFFICIENCYCALCULATOR_RMSY_RANGE5_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE5_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE5_MAX = 7.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE6 = "calibration-cdl-apr2014-Mn-Cr-12kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE6_LOW = 4.9
const EFFICIENCYCALCULATOR_ENERGY_RANGE6_HIGH = 6.9
const EFFICIENCYCALCULATOR_CHARGE_RANGE6_MIN = 3.5e5
const EFFICIENCYCALCULATOR_CHARGE_RANGE6_MAX = 6.0e5
const EFFICIENCYCALCULATOR_RMSY_RANGE6_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE6_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE6_MAX = 7.0
const EFFICIENCYCALCULATOR_ROOTTREE_RANGE7 = "calibration-cdl-apr2014-Cu-Ni-15kV"
const EFFICIENCYCALCULATOR_ENERGY_RANGE7_LOW = 6.9
const EFFICIENCYCALCULATOR_ENERGY_RANGE7_HIGH = 10.0
const EFFICIENCYCALCULATOR_CHARGE_RANGE7_MIN = 5.9e5
const EFFICIENCYCALCULATOR_CHARGE_RANGE7_MAX = 1.0e6 
const EFFICIENCYCALCULATOR_RMSY_RANGE7_MIN = 0.0
const EFFICIENCYCALCULATOR_RMSY_RANGE7_MAX = 1.1
const EFFICIENCYCALCULATOR_LENGTH_RANGE7_MAX = 7.0

var e0 = vec2(0.0)
e0[0] = EFFICIENCYCALCULATOR_ENERGY_RANGE0_LOW
e0[1] = EFFICIENCYCALCULATOR_ENERGY_RANGE0_HIGH
var e1 = vec2(0.0)
e1[0] = EFFICIENCYCALCULATOR_ENERGY_RANGE1_LOW
e1[1] = EFFICIENCYCALCULATOR_ENERGY_RANGE1_HIGH
var energyRanges = @[e0,e1]#,e2,e3,e4,e5,e6,e7]
echo energyRanges
echo energyRanges[0][0]

var c0 = vec2(0.0)
c0[0] = EFFICIENCYCALCULATOR_CHARGE_RANGE0_MIN
c0[1] = EFFICIENCYCALCULATOR_CHARGE_RANGE0_MAX
var c1 = vec2(0.0)
c1[0] = EFFICIENCYCALCULATOR_CHARGE_RANGE1_MIN
c1[1] = EFFICIENCYCALCULATOR_CHARGE_RANGE1_MAX
var chargeRanges = @[c0,c1]#,c2,c3,c4,c5,c6,c7]

proc findEfficiencySetting(chipRegion : string, softwareEfficiency: float64, efficiencySetting : float64) : int = 

  # function to find the correct efficiency setting for each energy region and write
  # them to efficiencySetting array
  # chipRegion:         considered region on the chip (gold, silver, bronze...)
  # softwareEfficiency: the desired software efficiency from which the cut value
  #                     is deduced
  # efficiencySetting:  array which stores the cut values, which are obtained

  var success : int
  success = 0
  var i : int
  for i in countup(0, EFFICIENCYCALCULATOR_NUMBER_RANGES):
    echo i
    #success = success + generateLikelihoodMarlinDistribution(rootTrees[i],chipRegion,chargeRanges[i][0],chargeRanges[i][1],lengthMaxs[i],rmsYRanges[i][0],rmsYRanges[i][1])

    #efficiencySetting[i] = findLikelihoodMarlinCutValue(softwareEfficiency)
  return success


proc generateLikelihoodMarlinDistribution(rootTree, chipRegion : string, totalChargeMin : float64, totalChargeMax : float64, lengthMax : float64, rmsYMin : float64, rmsYMax : float64) : int =
  var 
    goldCutActive = false
    silverCutActive = false
    bronzeCutActive = false
    wholeChipActive = false

  case chipRegion
  of "region: gold":
    goldCutActive = true
  of "region: silver":
    silverCutActive = true
  of "region: bronze":
    bronzeCutActive = true
  of "region: goldPlusSilver":
    goldCutActive = true
    silverCutActive = true
  of "region: goldPlusSilverPlusBronze":
    goldCutActive = true
    silverCutActive = true
    bronzeCutActive = true
  of "region: chip":
    wholeChipActive = true
  else: echo "Error: Unknown chip region!"
  return 1

  #// run over all events and for each check all cuts and if they are met,
  #// add likelihood value of this event to the _likelihoodMarlinDistribution histogram
  #[var iEvent : int
  for iEvent in 0..]#

const FILE = "calibration-cdl.h5"
proc getSize(h5F : string) : seq[int] =
  var 
    file_id: hid_t
    status: herr_t
    size: ptr csize
    #size: hid_t
  file_id = H5Fopen(h5F,H5F_ACC_RDWR, H5P_DEFAULT)
  file_id = hid_t H5Pget_size(file_id, h5F,size)
  status = H5Fclose(file_id)
  #status = H5Fclose(size)

echo getSize(FILE)
#
  #var status = H5Pget_size
  #return status



#findEfficiencySetting(chipRegionBackgroundAndDataChannelOne,softwareEfficiencyChannelOne,efficiencySettingChannelOne)


#type
#    vector3 = ref object of Vec3

## things changed## calibration-cdl.h5

# VT4 -> VT3
# XRT Focal length 1600.0 -> 1500.0
# RAYTRACER_RADIUS_PIPE_CB_VT3 = 33.6 #mm from drawing #30.0 #mm (original)
# RAYTRACER_LENGTH_PIPE_VT3_XRT = 264.7 #mm from picture #198.2 #mm (original)
# RAYTRACER_RADIUS_PIPE_VT3_XRT = 25.0 #mm from drawing #35.0 #m (original)
