import strscans, streams, strutils, math, os, tables, sequtils, strformat, hashes, polynumeric

import numericalnim, ggplotnim

import readSolarModel

type
  HeaderLine = enum
    H1, H2

  OpTableHeader = object
    case kind: HeaderLine
    of H1: density: int
    of H2: discard

  ElementKind = enum
    # uses proton number as value
    eH = 1
    eHe = 2
    eC = 6
    eN = 7
    eO = 8
    eNe = 10
    eNa = 11
    eMg = 12
    eAl = 13
    eSi = 14
    eP = 15
    eS = 16
    eCl = 17
    eAr = 18
    eK = 19
    eCa = 20
    eSc = 21
    eTi = 22
    eV = 23
    eCr = 24
    eMn = 25
    eFe = 26
    eCo = 27
    eNi = 28

  DensityOpacity = object
    ## a helper object to store the energy dependency of the opacity for a given
    ## density
    energies: seq[float]
    opacities: seq[float]
    # a cubic spline interpolation function to get any `energy` from the given
    # `energies` and `opacities`
    interp: CubicSpline[float]

  OpacityFile = object
    fname: string
    element: ElementKind
    temp: int
    densityTab: Table[int, DensityOpacity]

proc parseTableLine(energy, opacity: var float, line: string) {.inline.} =
  ## parses the energy and opacity float values from `line` into `energy`
  ## (1st col) and `opacity` (2nd col) using `scanf`
  if line.scanf("$s$f$s$f", energy, opacity):
    discard
  elif line.scanf("$s$f", opacity):
    discard
  else:
    raise newException(ValueError, "Parsing opacity table in line " & $line & " failed!")

proc parseTableHeader(line: string, hKind: HeaderLine): OpTableHeader =
  ## parses a line of the header of the monochromatic opacity file
  ## `hKind` is either the first or second line
  case hKind
  of H1:
    # do stuff for line 1, if we need something from here
    result = OpTableHeader(kind: H1)
    if line.scanf("$s$i", result.density):
      # NOTE: we use the density directly as `int` for simplicity reason!
      discard
    else: raise newException(ValueError, "Could not parse header line 1: " & $line)
  of H2:
    # do stuff for line 2, if we need something from here
    result = OpTableHeader(kind: H2)

proc parseOpacityFile(path: string): OpacityFile =
  ## we parse the monochromatic opacity file using strscans
  ## - first we drop the first line as the file header. Information in this?
  ## - then read table header (3 lines)
  ## - then num lines
  let ds = newFileStream(path)
  if ds.isNil:
    raise newException(IOError, "Could not open file " & $path)
  var buf = newString(200)
  var idx = 0
  let fname = path.extractFilename
  result = OpacityFile(fname: fname,
                       element: ElementKind(fname[2 .. 3].parseInt),
                       temp: parseInt(fname[5 .. ^1])) #pow(10.0, parseFloat(fname[5 .. ^1]) / 40.0))
  echo result
  var
    energy: float
    opac: float
    densityOpacity: DensityOpacity
  while not ds.atEnd:
    if idx == 0:
      # skip file header
      discard ds.readLine(buf)
      inc idx
      continue
    # read table header, 3 lines
    discard ds.readLine(buf)
    inc idx
    let h1 = parseTableHeader(buf, H1)
    discard ds.readLine(buf)
    let h2 = parseTableHeader(buf, H2)
    inc idx
    discard ds.readLine(buf)
    var tableCount = buf.strip.parseInt
    tableCount = if tableCount == 0: 10000 else: tableCount
    inc idx

    # now parse the table according to table count
    densityOpacity = DensityOpacity(energies: newSeq[float](tableCount),
                                    opacities: newSeq[float](tableCount))
    for j in 0 ..< tableCount:
      discard ds.readLine(buf)
      parseTableLine(energy, opac, buf)
      if tableCount == 10000:
        # set energy manually to `j`, since we simply have 1 eV steps
        energy = float j + 1
      densityOpacity.energies[j] = energy
      densityOpacity.opacities[j] = opac
      inc idx
    # finalize densityOpacity by creating spline and adding to result
    densityOpacity.interp = newCubicSpline(densityOpacity.energies,
                                           densityOpacity.opacities)
    result.densityTab[h1.density] = densityOpacity
    inc idx
  ds.close()

proc lagEval(n : int, x : float) : float = 
  if(n>0):
    result  = (2.0 * n.float - 1.0 - x) * lagEval(n-1,x) - (1.0 - (1.0 / n.float))* lagEval(n-2,x)
  elif(n==1):
    result =  1.0 - x
  elif(n==0):
    result =  1.0 
  else: result= 0.0

proc lagDeriv(m : int, x : float) : float = 
  if(m>0):
    result = lagDeriv(m-1,x) - lagEval(m-1,x)
  elif(m==0):
    result = 0.0 
  else : result = 0.0

proc quadWeight(x : float): float = 
  let N = 5
  result = 1.0 / (x * pow(lagDeriv(N,x),2) )

proc inner_integral(t : float, y: float) : float = 
  result = (1.0/2.0) * ( ((y * y) / (t * t + y * y)) + ln( t * t + y * y ) )


proc quadFunc(x : float, y : float, w : float): float = 
  var up_lim = sqrt(x+w) + sqrt(x);
  var lo_lim = sqrt(x+w) - sqrt(x);
  result = inner_integral(up_lim, y) - inner_integral(lo_lim, y);




proc F(w : float, y : float) : float =
  let N = 5
  var 
    res_coef : seq[float]
    lcoef = newSeqWith(N + 1, newSeq[float](N + 1))
    weights_vec : seq[float]
  lcoef[0][0] = 1.0
  lcoef[1][0] = 1.0
  lcoef[1][1]  = -1.0 # coeffs of the first two polynomials
  for n in 2..N: # n-th polynomial
    lcoef[n][0] = 1.0  #constants of all Laguerre polynomials are = 1
    for i in 1..n: # i-th power of x in the n-th polynomial
      lcoef[n][i] = ( (2 * n - 1).float * lcoef[n-1][i] - lcoef[n-1][i-1] + (1 - n).float * lcoef[n-2][i] ) / n.float
  ## storing the coefficients of the N-th order polynomial only whose roots will be computed
  for i in 0 .. N:
    res_coef.add(lcoef[N][i]) 
  #echo res_coef
  let p = initPoly(res_coef)
  ## calculate roots (Nullstellen) of this as a vector
  let roots_vec = p.roots()
  #for n in 0..N:
  #  roots_vec.add(0.0)
    #echo roots(res_coef[n])

  var integral = 0.0

  for i in 0 ..< N:
    weights_vec.add(quadWeight(roots_vec[i]))
  

  for i in 0 .. <N:
    integral += weights_vec[i] * quadFunc(roots_vec[i], y, w)
  

  integral *= (1.0 / 2.0)
  
  result = integral
let opFile = parseOpacityFile(testF)

# let's check whether the calculation worked by plotting the opacity for this file
# using the interpolation function we create
var dfSpline: DataFrame
for d, op in pairs(opFile.densityTab):
  let xs = linspace(1.0, 10000.0, 1000)
  let ys = xs.mapIt(op.interp.eval(it))
  let df = seqsToDf({ "energy" : xs,
                      "opacity" : ys,
                      "density" : toSeq(0 ..< xs.len).mapIt(d) })
  if dfSpline.len == 0:
    dfSpline = df
  else:
    dfSpline.add df

# filter out all opacities > 1.0 so that we can see the lines
proc str(i: Value): Value = %~ $i
let dfFiltered = dfSpline.filter(f{"opacity" < 1.0}).mutate(f{"densityStr" ~ str("density")})
# and plot all interpolated density opacities

ggplot(dfFiltered, aes("energy", "opacity", color = "densityStr")) +
  geom_line() +
  legendPosition(x = 0.8, y = 0.0) +
  ggtitle(&"E / Opacity for T = {opFile.temp} K, element: {opFile.element}") +
  ggsave("energy_opacity_density.pdf")

# alternatively plot all data in a log y plot
#ggplot(dfSpline, aes("energy", "opacity", color = "density")) +
#  geom_line() +
#  ggtitle("Energy dependency of opacity at different densities") +
#  scale_y_log10() +
#  ggsave("energy_opacity_density_log.pdf")



## First lets access the solar model and calculate some necessary values
const solarModel = "./ReadSolarModel/resources/AGSS09_solar_model_stripped.dat"

var df = readSolarModel(solarModel)
df = df.filter(f{"Radius" <= 0.2})
echo df.pretty(precision = 10)

# to read a single column, e.g. radius:
#echo df1["Rho"].len

# now let's plot radius against temperature colored by density
ggplot(df, aes("Radius", "Temp", color = "Rho")) +
  geom_line() +
  ggtitle("Radius versus temperature of solar mode, colored by density") +
  ggsave("radius_temp_density.pdf")

var 
  n_Z = newSeqWith(df["Rho"].len, newSeq[float](29)) #29 elements
  n_e : float
  n_es : seq[int]
  n_eInt : int
  distNe : float
  distTemp : float
  temperature : int
  temperatures : seq[int]
const
  alpha = 1.0 / 137.0
  g_ae = 1e-13
  m_e_keV = 510.998 #keV
  e_charge = 1.0
  kB = 1.380649e-23
let atomicMass = [1.0078,4.0026,3.0160,12.0000,13.0033,14.0030,15.0001,15.9949,16.9991,17.9991,20.1797,22.9897,24.3055,26.9815,28.085,30.9737,32.0675,35.4515,39.8775,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.9380,55.845,58.9331,58.6934] #all the 29 elements from the solar model file
let elements = ["H1", "He4","He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
# var solarTable = readSolarModel(solarModel)
let amu = 1.6605e-24 #grams
echo df["Mass"][6]
for iRadius in 0..< df["Rho"].len:
  n_Z[iRadius][1] = (df[elements[0]][iRadius].toFloat / atomicMass[0]) * (df["Rho"][iRadius].toFloat / amu) # Hydrogen
  for iZmult in 1..3:
    if iZmult == 1:
      n_Z[iRadius][iZmult * 2] = ((df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat) / ((atomicMass[iZmult * 2 - 1] * df[elements[iZmult * 2 - 1]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)) + (atomicMass[iZmult * 2] * df[elements[iZmult * 2]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
    else: n_Z[iRadius][iZmult + 4] = ((df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat) / ((atomicMass[iZmult * 2 - 1] * df[elements[iZmult * 2 - 1]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)) + (atomicMass[iZmult * 2] * df[elements[iZmult * 2]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu) 
  n_Z[iRadius][8] = ((df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat) / ((atomicMass[7] * df[elements[7]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[8] * df[elements[8]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[9] * df[elements[9]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
  for iZ in 10..<29:
    n_Z[iRadius][iZ] = (df[elements[iZ]][iRadius].toFloat / atomicMass[iZ]) * (df["Rho"][iRadius].toFloat / amu) # The rest
  n_e = (df["Rho"][iRadius].toFloat/amu) * (1 + df[elements[0]][iRadius].toFloat/2) # (g/cm³ /g) = 1/cm³
  #echo log(parseFloat(solarTable["Temp"][iRadius]), 10.0) / 0.025
  for iTemp in 0..90:
    distTemp = log(df["Temp"][iRadius].toFloat, 10.0) / 0.025 - float(140 + 2 * iTemp)
    if abs(distTemp) <= 1.0: 
      temperature = 140 + 2 * iTemp
      #echo distTemp
  temperatures.add(temperature)
  for iNe in 0..17:
    distNe = log(n_e, 10.0) / 0.25 - float(74 + iNe * 2)
    if abs(distNe) <= 1.0: 
      n_eInt = 74 + iNe * 2
  n_es.add(n_eInt)
  echo df["Temp"][iRadius].toFloat
  echo n_Z[iRadius][26]
#echo temperatures



#echo n_es



proc hash(x: ElementKind): Hash = 
  var h: Hash = 0
  result = h !& int(x)
  result = !$result

var densities: HashSet[int]
# var opElements: array[ElementKind, seq[OpacityFile]]
var opElements: Table[ElementKind, Table[int, OpacityFile]]
# allows for: opElements[Z][temp].densityTab[n_e].interp(E)
for temp in toSet(temperatures):
  for Z in ElementKind:
    let testF = &"./OPCD_3.3/mono/fm{int(Z):02}.{temp}"
    #echo Z
    #echo existsFile(testF)
    if existsFile(testF):
      let opFile = parseOpacityFile(testF)
      for k in keys(opFile.densityTab):
        densities.incl k
      if Z notin opElements:
        opElements[Z] = initTable[int, OpacityFile]()
      opElements[Z][temp] = opFile
      #[for iE in energies:
        for R in 0..< df["Rho"].len:
          if temp == temperatures[R]:
            n_eInt = n_es[R]  
            var opacity = opFile.densityTab[n_eInt].interp.eval(iE)]#


#echo densities

let energies = linspace(1.0, 10000.0, 1112)

var absCoefs = newSeqWith(df["Rho"].len, newSeq[float](1112)) #29 elements
var emratesS = newSeqWith(df["Rho"].len, newSeq[float](1112))
var ironOp = newSeqWith(df["Rho"].len, newSeq[float](1112))
var ironOpE : seq[float]

echo df["Rho"].len
let noElement = @[3, 4, 5, 9, 15, 17, 19, 21, 22, 23, 27]
for R in 0..<df["Rho"].len:
  n_eInt = n_es[R]
  temperature = temperatures[R]
  for iE in energies:
    var sum = 0.0
    var absCoef = 0.0
    var n_e_keV = pow(10.0, (n_es[R].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion
    var temp_keV = pow(10.0, (temperatures[R].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
    var temp_K = pow(10.0, (temperatures[R].toFloat * 0.025))
    var energy_keV = iE * 0.001
    var energy_J = 1.60218e-16 * energy_keV
    var w = energy_keV / temp_keV #(energy_J / ( kB * temp_K)) * 500.0
    var iEindex = find(energies, iE) # slows down the code: change
    #if w < 1.0:
      #echo "w" 
      #echo w
    #echo (energy_keV / temp_keV)
    for Z in ElementKind:
      if int(Z) in noElement: #Phosphorus and some other elements also don't exist in opacity files Z=15, etc.
        continue
      var opacity = opElements[Z][temperature].densityTab[n_eInt].interp.eval(iE)
      var opacity_keV = opacity * 0.528e-10 * 0.528e-10 * 5.076142e9 * 5.076142e9 # correct conversion
      # opacities in atomic unit for lenth squared: 0.528 x10-8cm * 0.528 x10-8cm = a0² # 1 m = 1/1.239841336215e-9 1/keV and a0 = 0.528 x10-10m
      sum += opacity_keV * n_Z[R][int(Z)] * 7.683e-24 #ToDo: check if this is the right way to transform the atomic number density into keV
      
    #echo iE.toInt    
    var n_e_keV = pow(10.0, (n_es[R].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion
    var temp_keV = pow(10.0, (temperatures[R].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
    var iEindex = find(energies, iE)
    absCoef[R][iEindex] = sum * (exp(iE * 0.001 / temp_keV) - 1.0) # is in keV
    #iE is in eV 
    # if want to have absorbtion coefficient of a radius and energy: R = (r (in % of sunR) - 0.0015) / 0.0005
    # energy = energies[iEindex]
    
    ## Now it's left to calculate the emission rates and for that the compton emission rate will be calculated first
    var energy_keV = iE * 0.001
    var compton_emrate = (alpha * g_ae * g_ae * energy_keV * energy_keV * n_e_keV) / (3.0 * m_e_keV * m_e_keV * (exp(energy_keV / temp_keV) - 1.0)) # (keV³ / keV²) = keV
    ## And the bremsstrahlung emission rate
    var debye_scale = sqrt( (4.0 * PI * alpha / temp_keV) * (n_e_keV + n_Z[R][1] * 7.645e-24 + 4.0 * n_Z[R][2] * 7.645e-24 )) ## making the same approximation as for n_e calculation 
    var y = debye_scale / (sqrt( 2.0 * m_e_keV * temp_keV))
    var brems_emrate = (alpha * alpha * g_ae * g_ae * (4.0/3.0) * sqrt(PI) * n_e_keV * n_e_keV * exp(-(energy_keV / temp_keV)) * gauss.F((energy_keV / temp_keV), sqrt(2) * y)) / (sqrt(temp_keV) * pow(m_e_keV,3.5) * energy_keV)

#echo absCoef






when false:



  proc getOpacity(opH: seq[OpacityFile], T, n_e, E: float): float =
    # angenommen T ist die Form die in OpacityFile enthalten
    let opF = opH.filterIt(it.temp == T)
    let dOp = opF.denstiyTab[n_e]
    result = dOp.interp(E)

  