import glm/vec
import math
import random
# kompilieren und ausführen: nim cpp -r aEL.nim, nim c -r --threads:on --showAllMismatches:on aEL.nim # nim cpp -r --gc:boehm --verbosity:3 aEL.nim ##hdfview likelihood_2018_2.h5  

import plotly
import random
import sequtils, os, strutils
import nimhdf5
#import ingrid/[tos_helpers, likelihood, ingrid_types]






##################rayTracer###############################


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
  
  var misalignmentDetector = vec3(0.0)
  misalignmentDetector[0] = detectorMisalignmentX
  misalignmentDetector[1] = detectorMisalignmentY
  misalignmentDetector[2] = 0

  var centerSun = vec3(0.0)
  centerSun[0] = 0
  centerSun[1] = 0
  centerSun[2] = RAYTRACER_DISTANCE_SUN_EARTH
  
  centerSun = centerSun + misalignmentSun

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

  var 
    radiusPipeVT3XRT = RAYTRACER_RADIUS_PIPE_VT3_XRT
    distanceCBAxisXRTAxis = RAYTRACER_DISTANCE_AXIS_CB_AXIS_XRT


  var
    integralNormalisation = 0.0
    integralTotal = 0.0
    integralDetector = 0.0
    integralBronze = 0.0
    integralSilver = 0.0
    integralGold = 0.0

  var misalignment = (sunMisalignmentH != 0.0) or (sunMisalignmentV != 0.0) or (detectorMisalignmentX != 0.0) or (detectorMisalignmentY != 0.0)
  # am Anfang fehlt noch etwas
  var plotName = ""
  plotname = plotname & axionRadiationCharacteristic #hier müssen wir wieder was überlegen, da das Objekt in den Klammern nochnichtnachgebaut ist
  plotname = plotname & "-txrtXX-" 
  var xrtTransmission_path_int = xrtTransmissionAt10Arcmin * 100.0 + 0.5 #//das ist eine übergebene Variable -> ok
  #echo xrtTransmission_path_int
  plotname = plotname & $xrtTransmission_path_int

  if misalignment == true:
    plotName = plotname & "-" & $(sunMisalignmentH) & "'"
    plotName = plotname & "-" & $(sunMisalignmentV) & "'"
    plotName = plotname & "-" & $(detectorMisalignmentX) & "mm"
    plotName = plotname & "-" & $(detectorMisalignmentY) & "mm"

  plotName = plotname & "-axion-image"


  ######### some functions we're gonna need later##########
  # Now let's define a function to get us a random point on a disk, where we later can put the center of the
  # detector as center and the radius of the detector as radius
  # , vecdata_x : array[0..1000, float], vecdata_y : array[0..1000,float]

  proc prepareheatmap( numberofrows : int, numberofcolumns : int, start_x : float, stop_x : float, start_y : float, stop_y : float, data_X :seq, data_Y : seq, weight1 : seq, norm : float64) : any =
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
      if coord_X > 0:
        if coord_Y > 0:
          if coord_X < float(numberofrows):
            if coord_Y < float(numberofcolumns):
              heatmaptable[int(coord_X)][int(coord_Y)] = heatmaptable[int(coord_X)][int(coord_Y)] + 1*weight1[i]/norm 
    result = heatmaptable

  proc getMaxVal(table : any, numberofrows : int) : float =
    var maxVals : seq[float]
    var maxVal : float64
    for i in 0 ..< numberofrows:
      maxVals.add(max(table[i]))
    maxVal = max(maxVals)
    result = maxVal

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
      y = @[float32(CHIPREGIONS_GOLD_Y_MIN * 280.0 / 14.0), float32(CHIPREGIONS_GOLD_Y_MIN * 280.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 280.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 280.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MIN * 280.0 / 14.0)]
      x = @[float32(CHIPREGIONS_GOLD_X_MIN * 280.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 280.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 280.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 280.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 280.0 / 14.0)]
    let
      d4 = Trace[float32](mode: PlotMode.LinesMarkers, `type`: PlotType.ScatterGL, ys: y, xs : x) 
      
                     
    let
      layout = Layout(title: diagramtitle, width: 800, height: 800,
                      xaxis: Axis(title: "x-axis [mm]"),#,  range: (0.0, 14.0)),
                      yaxis: Axis(title: "y-axis [mm]"), autosize: false)
      p = Plot[float32](layout: layout, traces: @[d, d4])
    echo p.save()
    p.show()



   
  proc getRandomPointOnDisk(center: Vec3, radius:float64) : Vec3 =
    #randomize()
    var 
      x = 0.0
      y = 0.0
      r = radius * sqrt(random(1.0))  #random(1.0)#Zufallszahl zwischen 0 und 1 (hoffe), die Wurzel nach C++-Vorlage(wtf)
    # _randomGEnerator -> Circle(x,y,r)  is done through the following
      angle = 280 * random(1.0) #random angle
    x = cos(angle) * r
    y = sin(angle) * r
    var vector = vec3(x,y,0.0)
    vector = vector + center
    result = vector 



  proc getRandomPointFromSolarModel(center : Vec3, radius : float64 ) : Vec3 =
    var 
      x = 0.0
      y = 0.0
      z = 0.0
      r = radius * 1e-1 * random(1.0) #e-1
    #in case of the standard axion radiation, we use 1/100 of the solar radius as 
    #the origin of axion radiation. In that region we assume homogeneous emission
      angle1 = 280 * random(1.0)
      angle2 = 180 * random(1.0)
    x = cos(angle1) * sin(angle2) * r
    y = sin(angle1) * sin(angle2) * r
    z = cos(angle2) * r
    var vector = vec3(x,y,z)
    vector = vector + center
    result = vector

  #echo getRandomPointFromSolarModel(centerSun,radiusSun) #Bedenke: der Ursprung dieser Koordinaten ist die Erde
  #echo getRandomPointFromSolarModel(centerSun,radiusSun)
  #echo getRandomPointFromSolarModel(centerSun,radiusSun)
  #echo random(1.0)  #seed is 0.8244128746481869

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


  #get the x and y values from the run-file to compare them to our model

  proc getXandY( h5file : string , dsetgrp1 :string, numFirstRun : int, numLastRun : int, chip : string, xOrY : string) : seq[float] =
    #[var h5f2 = H5File( h5file, "r")
    echo h5f2
    var g1_name = "likelihood"
    var g1 = h5f2[g1_name.grp_str] 
    echo g1
    var dsetgrpb = g1]#
    #var dsetgrp : H5Group
    #for dsetgrpb in items(h5f2, start_path = "likelihood"):
    
    var
      grp_name : string
      run_name : string
      valuesX : seq[float]
      valuesY : seq[float]
    for i in 0 .. (numLastRun-numFirstRun):
      run_name = "/run_" & $(numFirstRun+i)
      grp_name = dsetgrp1 & run_name & "/" & chip#& "/" & $(grp) & "/" & chip# 

      #if contains(dsetgrpb, dsetgrp1 & run_name): #contains(h5f, "likelihood")
      try:  
        withH5( h5file, "r"):
              # open h5 file using template
          var
            energy = h5f[(grp_name / "energyFromCharge"), float64]
              #logL = h5file[(grp_name / "likelihood"), float32]
            centerX = h5f[(grp_name / "centerX"), float64]
            centerY = h5f[(grp_name / "centerY"), float64]
              #ecc = h5file[(grp_name / "Excentricity"), float32]
              #length = h5file[(grp_name / "Length"), float32]
              #charge = h5file[(grp_name / "TotalCharge"), float32]
              #rmsTrans = h5file[(grp_name / "RmsTransverse"), float32]
              #npix = h5file[(grp_name / "NumberOfPixels"), float32]
          for i in 0 .. energy.high:
            valuesX.add(centerX[i])
            valuesY.add(centerY[i])
      except :
        echo "KeyError!"
        #else: echo @[0.0]
    if xOrY == "X":
      result = valuesX
    elif xOrY == "Y": 
      result = valuesY
    else: return @[0.0]

  ############done with the functions, let's use them############
  
  var pointdataX : seq[float] 
  var pointdataY : seq[float]
  var weights : seq[float]
  
  
  for iExitCB in 1..numberOfPointsEndOfCB:
    #echo iExitCB
    var pointExitCBMagneticField = vec3(0.0)
    pointExitCBMagneticField = getRandomPointOnDisk(centerExitCBMagneticField, radiusCB)
    for iSun in 1..numberOfPointsSun:
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
      var intersectsCB = false
        
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
      var pointExitPipeVT3XRT = vec3(0.0)
        
      if (not lineIntersectsCircle(pointExitCB,pointExitPipeCBVT3,centerExitPipeVT3XRT,radiusPipeVT3XRT,pointExitPipeVT3XRT)) : continue

      pointExitPipeVT3XRT = pointExitCB + ((centerExitPipeVT3XRT[2] - pointExitCB[2]) / (pointExitPipeCBVT3 - pointExitCB)[2]) * (pointExitPipeCBVT3 - pointExitCB)
        
      #pointExitPipeVT3XRT = pointInSun - pointExitCBMagneticField 
      #pointExitPipeVT3XRT = pointExitPipeVT3XRT/ sqrt(pointExitPipeVT3XRT[0]*pointExitPipeVT3XRT[0] + pointExitPipeVT3XRT[1]*pointExitPipeVT3XRT[1] + pointExitPipeVT3XRT[2]*pointExitPipeVT3XRT[2])
      #pointExitPipeVT3XRT = pointExitCBMagneticField + pointExitPipeVT3XRT * (RAYTRACER_LENGTH_COLDBORE + RAYTRACER_LENGTH_PIPE_CB_VT3 + RAYTRACER_LENGTH_PIPE_VT3_XRT)

      var vectorBeforeXRT = vec3(0.0)
      vectorBeforeXRT = pointExitPipeVT3XRT - pointExitCB

        
      ###################von CB (coldbore(pipe in Magnet)) zum XRT (XrayTelescope)#######################

      var pointEntranceXRT = vec3(0.0)
      pointEntranceXRT[0] = pointExitPipeVT3XRT[0]
      pointEntranceXRT[1] = pointExitPipeVT3XRT[1] + distanceCBAxisXRTAxis
      pointEntranceXRT[2] = pointExitPipeVT3XRT[2]

      var 
        angle = (arccos(vectorBeforeXRT[2]/sqrt(vectorBeforeXRT[0]*vectorBeforeXRT[0]+vectorBeforeXRT[1]*vectorBeforeXRT[1]+vectorBeforeXRT[2]*vectorBeforeXRT[2])))  # Here we want to adress theta, the polar angle, which should be the second entrance of the vector
        r_x = pointEntranceXRT[0]
        r_y = pointEntranceXRT[1]
        theta_x = vectorBeforeXRT[0] / vectorBeforeXRT[2]
        theta_y = vectorBeforeXRT[1] / vectorBeforeXRT[2]
        theta_x_prime = theta_x - ( r_x / RAYTRACER_FOCAL_LENGTH_XRT) 
        theta_y_prime = theta_y - ( r_y / RAYTRACER_FOCAL_LENGTH_XRT)
        
      var vectorAfterXRT = vec3(0.0)
      vectorAfterXRT[0] = sin(theta_x_prime) * 100.0
      vectorAfterXRT[1] = sin(theta_y_prime) * 100.0
      vectorAfterXRT[2] = 100.0
        
      var vectorAfterXRTPolar = vec3(0.0)  #(r,theta,phi)
      vectorAfterXRTPolar[0] = sqrt(vectorAfterXRT[0]*vectorAfterXRT[0]+vectorAfterXRT[1]*vectorAfterXRT[1]+vectorAfterXRT[2]*vectorAfterXRT[2])
      vectorAfterXRTPolar[1] = radToDeg(arccos(vectorAfterXRT[0]/vectorAfterXRTPolar[0]))
      vectorAfterXRTPolar[2] = radToDeg(arctan2(vectorAfterXRT[2],vectorAfterXRT[1]))

      vectorAfterXRTPolar[1] = 90.0 - vectorAfterXRTPolar[1] #this is the pitch angle
      var p = vectorAfterXRTPolar[1]
      vectorAfterXRTPolar[2] = 90.0 - vectorAfterXRTPolar[2] #this is the yaw angle
      var ya = vectorAfterXRTPolar[2]
      var transmissionTelescopePitch = 0.0008*p*p*p*p + 1e-04*p*p*p - 0.4489*p*p - 0.3116*p + 96.787
      #echo "now" 
      #echo angle
      #echo vectorAfterXRTPolar


      var centerDetectorWindow = vec3(0.0)
      centerDetectorWindow[0] = 0.0
      centerDetectorWindow[1] = 0.0
      centerDetectorWindow[2] = pointEntranceXRT[2] + RAYTRACER_FOCAL_LENGTH_XRT + RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW

      var lambda_0 = ( centerDetectorWindow[2] - pointEntranceXRT[2] ) / vectorAfterXRT[2]
      var pointDetectorWindow = vec3(0.0)
      pointDetectorWindow = pointEntranceXRT + lambda_0 * vectorAfterXRT
      pointDetectorWindow = pointDetectorWindow - misalignmentDetector

      
      #var weight = (transmissionTelescopePitch* (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) ) 
      var weight = ( telescopeTransmission(angle,xrtTransmissionAt10Arcmin) * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) )
      integralTotal = integralTotal + weight
        
      ###detector COS has (0/0) at the bottom left corner of the chip

      pointDetectorWindow[0] = pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X
      pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y
        
      pointdataX.add(pointDetectorWindow[0])
      pointdataY.add(pointDetectorWindow[1])
      weights.add(weight)


      #echo pointDetectorWindow

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
        
      #if(withinWindow){  image->Fill(pointDetectorWindow[0],pointDetectorWindow[1],weight)
  #echo pointdataX[0]
  # get the heatmaps out of the sequences of data X and data Y, first for the amount of data in one pixel
  # compared to the overall amount and then the data in one pixel compared to the maximal amount of data in any pixel   
  var 
    beginX = 0.0 # distanceCBAxisXRTAxis * 0.01
    endX = 14.0 #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 + distanceCBAxisXRTAxis * 0.01
    endY = 14.0 + distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareheatmap(280,280,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,integralNormalisation)#colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareheatmap(280,280,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,1.0)
  #echo heatmaptable2 #= 5417.0
  echo getMaxVal(heatmaptable2, 280)
  var heatmaptable3 = prepareheatmap(280,280,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,getMaxVal(heatmaptable2, 280)) # if change number of rows: has to be in the maxVal as well
  echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
  echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]

  #echo drawfancydiagrams("AxionModelFluxfraction", heatmaptable1, 280) 
  echo drawfancydiagrams("AxionModelProbability", heatmaptable3, 280) #Probabilities, that a photon, that hits a certain pixel could originate from an Axion, if the highest is 100%
  echo integralNormalisation # number of hits before the setup
  echo pointdataX.len # number of hits after the setup

  # get the heatmap of the data of a run for comparison 


  var FILE = "likelihood_2018_2.h5"
  var dataValuesX = getXandY( FILE,"likelihood",242, 306, "chip_3","X")
  var dataValuesY = getXandY( FILE,"likelihood",242, 306,"chip_3","Y")

  var weightData : seq[float]
  var weightProb : seq[float]
  echo dataValuesX.len
  for i in 0 ..< dataValuesX.len:
    weightData.add(1.0) 
    var X = int(dataValuesX[i]*20.0)
    var Y = int(dataValuesY[i]* 20.0)
    #echo heatmaptable3[X][Y] * 10.0
    weightProb.add(heatmaptable3[X][Y] * 100.0)
  #echo dataValuesX
  #echo dataValuesY
  var heatmaptable4 = prepareheatmap(280,280,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightData,1.0)
  #echo drawfancydiagrams("AxionModelDataRaw", heatmaptable4, 280) # the normal data of a run of chip 3 
  var heatmaptable5 = prepareheatmap(280,280,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightProb,1.0)
  #echo drawfancydiagrams("AxionModelProbability in %", heatmaptable5, 280) #the probability distribution of being an axion of the data of a run of chip 3

  

  #### now let's get the fluxfraction of the gold region by getting the weight of each event (probability of transition, dependend on the XRay-telescope transmission and the lenth of the 
  # path the particle would have traveled through the magnet) and divide it through the number of all events

  
  fluxFractionTotal    = integralTotal    / integralNormalisation
  fluxFractionDetector = integralDetector / integralNormalisation
  fluxFractionBronze   = integralBronze   / integralNormalisation
  fluxFractionSilver   = integralSilver   / integralNormalisation
  fluxFractionGold     = integralGold     / integralNormalisation

  echo "Flux fraction for the gold region:"
  echo getFluxFraction("region: gold")





var radiationCharacteristic : string ##axionRadiation::characteristic radiationCharacteristic(axionRadiation::characteristic::sar);
radiationCharacteristic = "axionRadiation::characteristic::sar"
var coldboreBlockedLength : float64 
coldboreBlockedLength = 0.0
var detectorWindowAperture : float64 
detectorWindowAperture = 14.0 #mm

echo calculateFluxFractions(radiationCharacteristic, xrtTransmissionAt10Arcmin, detectorWindowAperture, 0.0, 0.0, 0.0, 0.0, coldboreBlockedLength) # radiationCharacteristic = "axionRadiation::characteristic::sar"

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
