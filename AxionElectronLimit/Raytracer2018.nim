import glm/vec
import math
import random
# kompilieren und ausführen: nim cpp -r aEL.nim, nim c -r --threads:on --showAllMismatches:on aEL.nim # nim cpp -r --gc:boehm --verbosity:3 aEL.nim ##hdfview likelihood_2018_2.h5  

import plotly
import random
import sequtils, os, strutils
import nimhdf5
import chroma
#import ingrid/[tos_helpers, likelihood, ingrid_types]






##################rayTracer###############################

# axion conversion probability is at 1 for whole magnet and perpendicular to magnetic field
# put in XRT effeciency
# mirrors of XRT
# do photons fly the exact same path as the axion before?
# energy of the axiond and telescope transmission efficiency for that
#degToRad(angle has to be done in Raytracer2014 for cos and sin



################################
# VARIABLES from rayTracer.h
const RAYTRACER_DISTANCE_SUN_EARTH =  1.5e14  #mm #ok
const RAYTRACER_RADIUS_SUN = 6.9e11 #mm #ok
const RAYTRACER_RADIUS_COLDBORE = 21.5 #mm #ok
const RAYTRACER_LENGTH_COLDBORE = 9756.0 #mm half B field to end of CB #ok
const RAYTRACER_LENGTH_COLDBORE_9T = 9260.0 #mm half B field to half B field #ok
const RAYTRACER_LENGTH_PIPE_CB_VT3 = 2571.5 #mm should stay the same #from beam pipe drawings #ok
const RAYTRACER_RADIUS_PIPE_CB_VT3 = 39.64 #30.0 #mm smallest aperture between end of CB and VT4 # 43 mm but only 85% of area: 39,64mm #ok
const RAYTRACER_LENGTH_PIPE_VT3_XRT = 150 #mm from drawings #198.2 #mm from XRT drawing #ok
const RAYTRACER_RADIUS_PIPE_VT3_XRT = 35.0#25.0 #mm from drawing #35.0 #m irrelevant, large enough to not loose anything # needs to be mm #ok
const RAYTRACER_FOCAL_LENGTH_XRT = 1485.0 #mm is from llnl XRT https://iopscience.iop.org/article/10.1088/1475-7516/2015/12/008/pdf #1600.0 #mm was the Telescope of 2014 (MPE XRT) also: Aperatur changed #ok
const RAYTRACER_DISTANCE_AXIS_CB_AXIS_XRT = 0.0#62.1#58.44 #mm from XRT drawing #there is no difference in the axis even though the picture gets transfered 62,1mm down, but in the detector center
const RAYTRACER_DISTANCE_FOCAL_PLANE_DETECTOR_WINDOW = 20.0 #mm #no change, because don't know # is actually -10.0 mm
const numberOfPointsEndOfCB = 4000
const numberOfPointsSun = 4000
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
  ## first functions for the graphs later, that store the data in heatmaps and then give them out with plotly


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
      if coord_X >= 0:
        if coord_Y >= 0:
          if coord_X <= float(numberofrows):
            if coord_Y <= float(numberofcolumns):
              heatmaptable[int(coord_X)][int(coord_Y)] = heatmaptable[int(coord_X)][int(coord_Y)] + 1*weight1[i]/norm 
    result = heatmaptable

  proc getMaxVal(table : any, numberofrows : int) : float =
    var maxVals : seq[float]
    var maxVal : float64
    for i in 0 ..< numberofrows:
      maxVals.add(max(table[i]))
    maxVal = max(maxVals)
    result = maxVal

  proc getLenght(table : any, numberofrows : int) : int =
    var lengths : seq[float]
    var length : int
    for i in 0 ..< numberofrows:
      lengths.add(max(table[i]))
    length = lengths.len
    result = length

  proc drawfancydiagrams(diagramtitle : string, objectstodraw : any, width : int) : float =
    let
      # The GL heatmap is also supported as HeatMapGL
      d = Trace[float32](mode: PlotMode.Lines, `type`: PlotType.HeatMap)
    
    
    # fill data for colormap with random values. The data needs to be supplied
    # as a nested seq.
    #
    #
    #
    #
    #
    #[var
      x2 : seq[float32]
      y2 : seq[float32]
    const color_choice = @[Color(r: 0.9, g: 0.1, b: 0.1, a: 1.0)]
    var 
      colors : seq[Color] #colors = newSeqWith(width, newSeq[Color](width))
      sizes : seq[float32] #sizes = newSeqWith(width, newSeq[float64](width))
      vals : seq[float32]

    for i in 0 ..< width:
      for j in 0 ..< width:
        if objectstodraw[i][j] == 0.0:
          colors.add(color_choice[0])#colors[i][j] = color_choice[0]
          sizes.add(20.0)#sizes[i][j] = 2.0
          vals.add(0.0)
          x2.add(float32(i))
          y2.add(float32(j))]#
        #else:
          #colors[i][j] = ColorMap.Viridis
    d.colormap = ColorMap.Viridis    
    d.zs = newSeqWith(width, newSeq[float32](width))
    for x in 0 ..< width:
      for y in 0 ..< width:
        if x < width:
          if y < width:
            if x > 0:
              if y > 0:
                d.zs[x][y] = objectstodraw[x][y] 
    #[let d2 = Trace[float32](mode: PlotMode.Markers, `type`: PlotType.HeatMap,ys: y2, xs : x2 ) 
    d2.marker = Marker[float32](size: sizes, color: colors)]##d2.marker = Marker[float32](size: sizes, colorVals: vals, colorMap: ColorMap.Jet)#

    const 
      y = @[float32(CHIPREGIONS_GOLD_Y_MIN * 3500.0 / 14.0), float32(CHIPREGIONS_GOLD_Y_MIN * 3500.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 3500.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MAX * 3500.0 / 14.0),float32(CHIPREGIONS_GOLD_Y_MIN * 3500.0 / 14.0)]
      x = @[float32(CHIPREGIONS_GOLD_X_MIN * 3500.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 3500.0 / 14.0),float32(CHIPREGIONS_GOLD_X_MAX * 3500.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 3500.0 / 14.0), float32(CHIPREGIONS_GOLD_X_MIN * 3500.0 / 14.0)]
    let
      d4 = Trace[float32](mode: PlotMode.LinesMarkers, `type`: PlotType.ScatterGL, ys: y, xs : x) 
    #d.marker = Marker[float64](size: sizes, Color(r: 0.9, g: 0.1, b: 0.1, a: 1.0))
      
                     
    let
      layout = Layout(title: diagramtitle, width: 800, height: 800,
                      xaxis: Axis(title: "x-axis [mm]"),#,  range: (0.0, 14.0)),
                      yaxis: Axis(title: "y-axis [mm]"), autosize: false)
      p = Plot[float32](layout: layout, traces: @[d, d4])
    echo p.save()
    p.show()

  # Now let's define a function to get us a random point on a disk, where we later can put the center of the
  # detector as center and the radius of the detector as radius

   
  proc getRandomPointOnDisk(center: Vec3, radius:float64) : Vec3 =
    #randomize()
    var 
      x = 0.0
      y = 0.0
      r = radius * sqrt(random(1.0))  #random(1.0)#why sqrt?
    # _randomGEnerator -> Circle(x,y,r)  is done through the following ### gives difference, since the root function Circle probably uses a differen random algorithm
      angle = 360 * random(1.0)  #random angle
    x = cos(degToRad(angle)) * r
    y = sin(degToRad(angle)) * r
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
      angle1 = 360 * random(1.0)
      angle2 = 180 * random(1.0)
    x = cos(degToRad(angle1)) * sin(degToRad(angle2)) * r
    y = sin(degToRad(angle1)) * sin(degToRad(angle2)) * r
    z = cos(degToRad(angle2)) * r
    var vector = vec3(x,y,z)
    vector = vector + center
    result = vector

  #echo getRandomPointFromSolarModel(centerSun,radiusSun) #Bedenke: der Ursprung dieser Koordinaten ist die Erde
  #echo getRandomPointFromSolarModel(centerSun,radiusSun)
  #echo getRandomPointFromSolarModel(centerSun,radiusSun)
  #echo random(1.0)  #seed is 0.8244128746481869

  # Now a function to see if the linesfrom the sun will actually intersect the circle area from the magnet entrance (called coldbore) etc.

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

  # Also a function to know if the line intersected at least the whole magnet, and then only once, because else the axions would have just flown through

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

  # telescope transmission of the old XRT, is not used in 2018

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

  #A function to create the coordinates of the mirrors of the telescope (they're made of glass) to substract from the overall picture, because Xrays in that range don't pass through glass

  proc circleEdges( theR1 : seq[float]): any =
    const 
      sizeViewfield = 45.5 #mm
      d = 83.0 #mm
      
    var 
      circleEdge = vec3(0.0)
      circleEdges = newSeqWith(1401, newSeq[float64](1401))#: seq[seq[float]] #
      #circleEdgesY : seq[float]
      #circleEdgesX : seq[float]

    for R1 in theR1:
      for phi in 16500 ..< 19500:
        for dgl in 0 .. 20:
          circleEdge[0]= (R1 + (float(dgl)*0.01)) * sin(degToRad(float(phi)*0.01)) 
          circleEdge[1]= (R1 + (float(dgl)*0.01)) * cos(degToRad(float(phi)*0.01)) + d
          var coord_X = floor(circleEdge[0] / (sizeViewfield/1400.0)) + 700
          var coord_Y = floor(circleEdge[1] / (sizeViewfield/1400.0)) + 700
          if coord_X >= 0.0:
            if coord_Y >= 0.0:
              if coord_X <= 1400.0:
                if coord_Y <= 1400.0:
                    circleEdges[int(coord_X)][int(coord_Y)] = circleEdges[int(coord_X)][int(coord_Y)] + 1
          #circleEdgesX.add(circleEdge[0])
          #circleEdgesY.add(circleEdge[1])
    #if xOrY == "X":
        #result = circleEdgesX
    #elif xOrY == "Y": 
        #result = circleEdgesY
    result = circleEdges 

  proc getPixelValue(intersects : Vec3): Vec3 =
    const sizeViewfield = 45.5 #mm
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
    var r_xy_intersect = sqrt(intersect[0] * intersect[0] + (intersect[1]+d) * (intersect[1]+d))
    if r_xy_intersect <= R1 and r_xy_intersect >= (prevR1) :
      return true
    else: 
      return false

  #get the x and y values from the run-file to compare them to our model

  proc getXandY( h5file : string , dsetgrp1 :string, numFirstRun : int, numLastRun : int, chip : string, xOrY : string) : seq[float] =

    
    var
      grp_name : string
      run_name : string
      valuesX : seq[float]
      valuesY : seq[float]
    for i in 0 .. (numLastRun-numFirstRun):
      run_name = "/run_" & $(numFirstRun+i)
      grp_name = dsetgrp1 & run_name & "/" & chip#& "/" & $(grp) & "/" & chip# 

      # if the run file with the run number does not exist, it is thrown away in this next step
      try:  
        withH5( h5file, "r"):
              # open h5 file using template
          var
            energy = h5f[(grp_name / "energyFromCharge"), float64]
            centerX = h5f[(grp_name / "centerX"), float64]
            centerY = h5f[(grp_name / "centerY"), float64]

          for i in 0 .. energy.high:
            if centerX[i] < 14.0:
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
  var pointdataXBefore : seq[float] 
  var pointdataYBefore : seq[float]
  var weights : seq[float]
  const 
    allR1 = @[60.7095, 63.006 , 65.606 , 68.305 , 71.105 , 74.011 , 77.027 , 80.157 , 83.405 , 86.775 , 90.272 , 93.902 , 97.668 , 101.576 , 105.632]
    circleX = circleEdges( allR1)
    circleY = circleEdges( allR1)
    circleTotal = circleEdges( allR1)
    d = 83.0 #mm
  var pixvalsX : seq[float] 
  var pixvalsY : seq[float]
  
  #echo min(circleX)
  #echo min(circleY)
  #echo circleTotal

  # count through all points in the magnet end (coldbore) and then in the sun (set to 1000 points), to connect them and get 1000000 lines (axion traces)

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
      pointEntranceXRT[0] = pointExitPipeVT3XRT[0] #- distanceCBAxisXRTAxis
      pointEntranceXRT[1] = pointExitPipeVT3XRT[1] 
      pointEntranceXRT[2] = pointExitPipeVT3XRT[2]
      #echo pointEntranceXRT
      #echo 
      
      if lineIntersectsCircleEdge(circleTotal, getPixelValue(pointEntranceXRT)): continue
        
      #echo getPixelValue(pointEntranceXRT)
      #echo lineIntersectsCircleEdge(circleTotal, getPixelValue(pointEntranceXRT))
      # x and y is interchanged from this point on since it somehow is interchanged in the picture
      if  pointEntranceXRT[0] <= 1.0 and pointEntranceXRT[0] >= -1.0: continue #there is a 2mm wide graphit block between each glass mirror, to seperate them in the middle of the x-Ray telescope
      
      var 
        vectorAfterXRTCircular = vec3(0.0)
        radius1 = sqrt(pointEntranceXRT[0] * pointEntranceXRT[0] + (pointEntranceXRT[1]+d) * (pointEntranceXRT[1]+d))
        phi_radius = arctan2(pointEntranceXRT[0],(pointEntranceXRT[1]+d)) #arccos((pointEntranceXRT[1]+d) / radius1)
        alpha = arctan(radius1 / RAYTRACER_FOCAL_LENGTH_XRT)
      vectorAfterXRTCircular[0] = radius1
      vectorAfterXRTCircular[1] = phi_radius
      vectorAfterXRTCircular[2] = alpha

      
      
      #for i in 1 .. < allR1.len:
        #echo lineIntersectsArea( (allR1[i] + 0.2), allR1[i], pointEntranceXRT)
        #if lineIntersectsArea( (allR1[i] + 0.2), allR1[i], pointEntranceXRT):continue
      

      var valuesPix = getPixelValue(pointEntranceXRT)
      pointdataXBefore.add(pointEntranceXRT[0])
      pointdataYBefore.add(pointEntranceXRT[1])
      pixvalsX.add(valuesPix[0])
      pixvalsY.add(valuesPix[1])
      #[var r_xy_XRT = sqrt(pointEntranceXRT[0] * pointEntranceXRT[0] + (pointEntranceXRT[1]+d) * (pointEntranceXRT[1]+d))
      #echo r_xy_XRT
      #echo lineIntersectsArea( allR1[2], 0.0, pointEntranceXRT)
      if lineIntersectsArea( allR1[1], 0.0, pointEntranceXRT):
        pointEntranceXRT[0] = sqrt((r_xy_XRT - allR1[0]) * (r_xy_XRT - allR1[0]) - pointEntranceXRT[1] * pointEntranceXRT[1])
        pointEntranceXRT[1] = sqrt((r_xy_XRT - allR1[0]) * (r_xy_XRT - allR1[0]) - pointEntranceXRT[0] * pointEntranceXRT[0])
      elif lineIntersectsArea( allR1[2], allR1[1], pointEntranceXRT):

        pointEntranceXRT[0] = sqrt((r_xy_XRT - (allR1[1] + 0.2)) * (r_xy_XRT - (allR1[1] + 0.2)) - pointEntranceXRT[1] * pointEntranceXRT[1])
        pointEntranceXRT[1] = sqrt((r_xy_XRT - (allR1[1] + 0.2)) * (r_xy_XRT - (allR1[1] + 0.2)) - pointEntranceXRT[0] * pointEntranceXRT[0])
      elif lineIntersectsArea( allR1[3], allR1[2], pointEntranceXRT):
        pointEntranceXRT[0] = sqrt((r_xy_XRT - (allR1[2] + 0.2)) * (r_xy_XRT - (allR1[2] + 0.2)) - pointEntranceXRT[1] * pointEntranceXRT[1])
        pointEntranceXRT[1] = sqrt((r_xy_XRT - (allR1[2] + 0.2)) * (r_xy_XRT - (allR1[2] + 0.2)) - pointEntranceXRT[0] * pointEntranceXRT[0])
      elif lineIntersectsArea( allR1[4], allR1[3], pointEntranceXRT):
        pointEntranceXRT[0] = sqrt((r_xy_XRT - (allR1[3] + 0.2)) * (r_xy_XRT - (allR1[3] + 0.2)) - pointEntranceXRT[1] * pointEntranceXRT[1])
        pointEntranceXRT[1] = sqrt((r_xy_XRT - (allR1[3] + 0.2)) * (r_xy_XRT - (allR1[3] + 0.2)) - pointEntranceXRT[0] * pointEntranceXRT[0])
      elif lineIntersectsArea( allR1[4], allR1[3], pointEntranceXRT):
        pointEntranceXRT[0] = sqrt((r_xy_XRT - (allR1[3] + 0.2)) * (r_xy_XRT - (allR1[3] + 0.2)) - pointEntranceXRT[1] * pointEntranceXRT[1])
        pointEntranceXRT[1] = sqrt((r_xy_XRT - (allR1[3] + 0.2)) * (r_xy_XRT - (allR1[3] + 0.2)) - pointEntranceXRT[0] * pointEntranceXRT[0])
      else: continue]#
      #echo pointEntranceXRT
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
      var pointDetectorWindow = vec3(0.0)
      pointDetectorWindow[0] = pointDetectorWindowCircle[0] * sin(pointDetectorWindowCircle[1])
      pointDetectorWindow[1] = pointDetectorWindowCircle[0] * cos(pointDetectorWindowCircle[1])
      pointDetectorWindow[2] = 0.0
      
      #echo vectorAfterXRTCircular
      #echo pointDetectorWindow

      #[var lambda_0 = ( centerDetectorWindow[2] - pointEntranceXRT[2] ) / vectorAfterXRT[2]
      
      pointDetectorWindow = pointEntranceXRT + lambda_0 * vectorAfterXRT
      pointDetectorWindow = pointDetectorWindow - misalignmentDetector]#
      
        
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
        transmissionMagnet = cos(ya) * 1.0 # this is the transformation probability of an axion into a photon, if an axion flying straight through the magnet had one of 100%, angular dependency of the primakoff effect



      var weight = (transmissionTelescopePitch*transmissionTelescopeYaw* transmissionMagnet * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) ) #transmission probabilities times time the axion spend in the magnet
      #var weight = ( telescopeTransmission(angle,xrtTransmissionAt10Arcmin) * (pathCB * pathCB / RAYTRACER_LENGTH_COLDBORE_9T / RAYTRACER_LENGTH_COLDBORE_9T) )
      integralTotal = integralTotal + weight
     
      ###detector COS has (0/0) at the bottom left corner of the chip

      pointDetectorWindow[0] = pointDetectorWindow[0] + CHIPREGIONS_CHIP_CENTER_X
      pointDetectorWindow[1] = pointDetectorWindow[1] + CHIPREGIONS_CHIP_CENTER_Y
      #echo pointDetectorWindow
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
    beginX = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endX = 14.0 #- distanceCBAxisXRTAxis * 0.01
    beginY = 0.0 #- distanceCBAxisXRTAxis * 0.01
    endY = 14.0 #- distanceCBAxisXRTAxis * 0.01
  var heatmaptable1 = prepareheatmap(3500,3500,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,integralNormalisation)#colour scale is now the number of points in one pixel divided by the the number of all events
  var heatmaptable2 = prepareheatmap(3500,3500,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,1.0)
  #echo heatmaptable2 #= 5417.0
  echo getMaxVal(heatmaptable2, 3500)
  var heatmaptable3 = prepareheatmap(3500,3500,beginX,endX,beginY,endY,pointdataX,pointdataY,weights,getMaxVal(heatmaptable2, 3500)) # if change number of rows: has to be in the maxVal as well
  echo "Probability of it originating from an axion if a photon hits at x = 5,3mm and y = 8,4mm (in this model):"
  echo (heatmaptable3[53][84]) * 100.0  #echo heatmaptable3[x][y]

  #echo drawfancydiagrams("AxionModelFluxfraction", heatmaptable1, 3500) 
  echo drawfancydiagrams("AxionModelProbability", heatmaptable3, 3500) #Probabilities, that a photon, that hits a certain pixel could originate from an Axion, if the highest is 100%
  echo integralNormalisation # number of hits before the setup
  echo pointdataX.len # number of hits after the setup

  # get the heatmap of the data of a run for comparison 
  echo "here"
  echo weights.len


  var FILE = "likelihood_2018_2_all.h5"
  var dataValuesX = getXandY( FILE,"likelihood",240, 306, "chip_3","X")
  var dataValuesY = getXandY( FILE,"likelihood",240, 306,"chip_3","Y")
  
  var weightData : seq[float]
  
  var weightProb : seq[float]
  echo getLenght(heatmaptable3, 3500)
  #echo dataValuesX.len
  for i in 0 ..< dataValuesX.len:
    weightData.add(1.0) 
    var X = int(floor(dataValuesX[i]*100.0))
    var Y = int(floor(dataValuesY[i]*100.0))
    #echo heatmaptable3[X][Y] * 10.0
    weightProb.add(heatmaptable3[X][Y] * 100.0)

  #echo dataValuesX
  #echo dataValuesY
  var heatmaptable4 = prepareheatmap(3500,3500,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightData,1.0)
  #echo drawfancydiagrams("AxionModelDataRaw", heatmaptable4, 3500) # the normal data of a run of chip 3 
  var heatmaptable5 = prepareheatmap(3500,3500,0.0,14.0,0.0,14.0,dataValuesX,dataValuesY,weightProb,1.0)
  #echo drawfancydiagrams("AxionModelProbability in %", heatmaptable5, 3500) #the probability distribution of being an axion of the data of a run of chip 3
  #[var weightData1 : seq[float]
  for i in 0 ..< circleX.len:
    weightData1.add(1.0)

  var heatmaptable6 = prepareheatmap(140,140,0.0,33.0,0.0,33.0,circleX,circleY,weightData1,1.0)
  echo drawfancydiagrams("Mirrors", heatmaptable6, 140)]#

  var weightData2 : seq[float]
  for i in 0 ..< pixvalsY.len:
    weightData2.add(1.0)

  var heatmaptable7 = prepareheatmap(1400,1400,0.0,1400.0,0.0,1400.0,pixvalsX,pixvalsY,weightData2,1.0)
  #echo drawfancydiagrams("Mirrors", heatmaptable7, 1400)
  #echo drawfancydiagrams("Mirrors",circleTotal , 1400)

  

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
