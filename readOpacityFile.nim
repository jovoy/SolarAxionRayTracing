import strscans, streams, strutils, math, os, tables, sequtils, strformat, hashes, polynumeric, macros

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

  OpacityFileKind = enum
    ofkOriginal, ofkNew

  OpacityFile = object
    fname: string
    element: ElementKind
    temp: int
    case kind: OpacityFileKind
    of ofkOriginal:
      densityTab: Table[int, DensityOpacity]
    of ofkNew:
      density: int
      densityOp: DensityOpacity

  ZTempDensity = tuple[Z: int, temp: int, density: int]

macro iterEnum(en: typed): untyped =
  result = nnkBracket.newTree()
  let impl = en.getImpl
  expectKind impl[2], nnkEnumTy
  var first = true
  for el in impl[2]:
    if first: 
      first = false
      continue
    result.add nnkPar.newTree(newLit(el[0].toStrLit.strVal), el[1])

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

proc convertLogUToE(logU, temp: float): float =
  ## converts a given energy in `logU` to a temperature in `eV`
  result = exp(logU) * pow(10.0, (temp * 0.025)) * 8.617e-8

proc parseDensityTab(ds: FileStream, temp: int, 
                     kind: OpacityFileKind,
                     tableCount = 1000): DensityOpacity =
  # now parse the table according to table count
  result = DensityOpacity(energies: newSeqOfCap[float](tableCount),
                          opacities: newSeqOfCap[float](tableCount))
  var 
    buf = newString(200)
    energy: float
    opac: float
    idx = 0
  while not ds.atEnd:
    discard ds.readLine(buf)
    parseTableLine(energy, opac, buf)
    case kind
    of ofkOriginal:
      if tableCount == 10000:
        # set energy manually to `j`, since we simply have 1 eV steps
        energy = float idx + 1
    of ofkNew:
      energy = convertLogUToE(energy, temp.float)
    result.energies.add energy
    result.opacities.add opac
    inc idx
    if kind == ofkOriginal and idx == tableCount:
      break
  # finalize densityOpacity by creating spline and adding to result
  result.interp = newCubicSpline(result.energies,
                                  result.opacities) 

proc parseOpacityNew(ds: FileStream, 
                         fname: string,
                         kind: static OpacityFileKind): OpacityFile =          
  var buf = newString(200)
  discard ds.readLine(buf)
  proc removeSuffix(s, suffix: string): string =
    result = s
    result.removeSuffix(suffix)
  let fnameSeq = fname.removeSuffix(".dat").split("_")
  let elStr = fnameSeq[2]
  let temp = parseInt(fnameSeq[3])
  result = OpacityFile(fname: fname,
                       kind: kind,
                       element: parseEnum[ElementKind]("e" & elStr),
                       temp: temp,
                       density: parseInt(fnameSeq[4]))
  result.densityOp = ds.parseDensityTab(temp = temp, kind = kind)

proc parseOpacityOriginal(ds: FileStream,
                          fname: string,
                          kind: static OpacityFileKind): OpacityFile =
  let temp = parseInt(fname[5 .. ^1])                        
  result = OpacityFile(fname: fname,
                       kind: kind,
                       element: ElementKind(fname[2 .. 3].parseInt),
                       temp: temp) #pow(10.0, parseFloat(fname[5 .. ^1]) / 40.0))
  var idx = 0
  var buf = newString(200)
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
    result.densityTab[h1.density] = ds.parseDensityTab(temp, kind, tableCount)
    inc idx                                             

proc parseOpacityFile(path: string, kind: OpacityFileKind): OpacityFile =
  ## we parse the monochromatic opacity file using strscans
  ## - first we drop the first line as the file header. Information in this?
  ## - then read table header (3 lines)
  ## - then num lines
  let ds = newFileStream(path)
  echo path
  if ds.isNil:
    raise newException(IOError, "Could not open file " & $path)
  let fname = path.extractFilename
  case kind
  of ofkOriginal:
    result = ds.parseOpacityOriginal(fname = fname,
                                     kind = ofkOriginal)
  of ofkNew:
    result = ds.parseOpacityNew(fname = fname,
                                    kind = ofkNew)
  ds.close()

proc readMeshFile(fname: string): seq[float] =
  var uValues : seq[float]
  var df = toDf(readCsv(fname))#, sep = ' ', header = "#"))
  for i in 0..10000:
    uValues.add(toFloat(df["u"][i]))
  result = uValues
  #[var
    f = open(fname) 
    line = ""
    ustring : string
    uValues : seq[float]

  while f.readLine(line):
    let ustring = line
    uValues.add(parseFloat(ustring))

  defer: f.close()
  return uValues]#

let lineNumbers = linspace(0.0, 10000.0, 10001)
const meshFile = "./OPCD_3.3/mono/fm01.mesh"
let dfMesh = readMeshFile(meshFile)

let spline = newCubicSpline(dfMesh, lineNumbers )

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

proc outer(x, w, y: float): float =
  let coeff = x * exp(-x * x)
  # wrap the inner call to have an `IntegrateFunc` kind
  #let fn = proc(x: float, optional: seq[float]): float =
    #result = inner(x, y)
  let
    frm = sqrt(x * x + w) - x
    to = sqrt(x * x + w) + x
  let integral = inner_integral(to, y) - inner_integral(frm, y)#adaptiveGauss(fn, frm, to) #check if analytical integral gives the same results: difference of e-15

  result = coeff * integral

proc fNew(w : float, y : float) : float =
  # integrate `outer` frm `0` to `inf`
  # rewrite via
  # int_a^infty dx f(x) = int_0^1 f(a + (1 - t) / t) / t^2
  # in our case a = 0
  # so express by wrapping `outer` in a new proc
  let fnToInt = proc(t: float, optional: seq[float]): float =
    #echo t
    if t != 0:
      result = outer(x = (1 - t) / t, w = w, y = y) / (t * t)
    else:
      # workaround singularity by adding epsilon
      result = outer(x = (1 - t) / (t + 1e-8), w = w, y = y) / (t * t)
  # and integrate that from 0 to 1
  result = adaptiveGauss(fnToInt, 0.0, 1.0)


proc quadFunc(x : float, y : float, w : float): float = 
  var up_lim = sqrt(x * x + w) + x # sqrt(x  +w) + sqrt(x)
  var lo_lim = sqrt(x * x + w) - x # sqrt(x  +w) - sqrt(x)
  result = inner_integral(up_lim, y) - inner_integral(lo_lim, y);




proc f(w : float, y : float) : float =
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
  var integral = 0.0
  for i in 0 ..< N:
    weights_vec.add(quadWeight(roots_vec[i]))

  for i in 0 ..< N:
    integral += weights_vec[i] * quadFunc(roots_vec[i], y, w)
  integral *= (1.0 / 2.0)
  result = integral


## The functions for all the parts of the emission rate ##
    
proc comptonEmrate(alpha : float, gae : float, energy : float, ne : float, me : float, temp : float) : float =
  result = (alpha * gae * gae * energy * energy * ne) / (3.0 * me * me * me * me * (exp(energy / temp) - 1.0))

proc bremsEmrate(alpha : float, gae : float, energy : float, ne : float, me : float, temp : float, w : float, y : float) : float =
  result = (alpha * alpha * gae * gae * 4.0 * sqrt(PI) * ne * ne * exp(- energy / temp) * fNew(w, sqrt(2.0) * y)) / (3.0 * sqrt(temp) * pow(me, 3.5) * energy)

proc term1(gae : float, energy : float, abscoef : float, echarge : float, me : float, temp : float) : float = 
  result = (gae * gae * energy * energy * abscoef) / (2.0 * echarge * echarge * me * me * (exp(energy / temp) - 1.0))

proc term2(alpha : float, gae : float, energy : float, ne : float, me : float, temp : float) : float =
  result = ((exp(energy / temp) - 2.0) * comptonEmrate(alpha, gae, energy, ne, me, temp)) / (2.0 * (exp(energy / temp) - 1.0))

proc freefreeEmrate(alpha : float, gae : float, energy : float, ne : float, me : float, temp : float, nzZ2 :float, w:float, y:float) : float =
  result = (fNew(w, y) * alpha * alpha * gae * gae  * 8.0 * sqrt(PI) * ne * nzZ2 * exp(-energy/temp)) / (3.0 * sqrt(2.0 * temp) * pow(me, 3.5) * energy)#*( ) * )


proc getFluxFraction(energies : seq[float], df : DataFrame, n_es : seq[int], temperatures : seq[int], emratesS : seq[seq[float]]) : seq[float] = 
  const
    alpha = 1.0 / 137.0
    g_ae = 1e-13 # Redondo 2013: 0.511e-10 
    m_e_keV = 510.998 #keV
    e_charge = sqrt(4.0 * PI * alpha)#1.0
    kB = 1.380649e-23
    r_sun = 6.957e11 #mm
    r_sunearth = 1.5e14 #mm
    hbar = 6.582119514e-25 # in GeV * s
    keV2cm = 1.97327e-8 # cm per keV^-1
    amu = 1.6605e-24 #grams
  var diff_fluxs : seq[float]
  let factor =   pow( r_sun * 0.1 / (keV2cm), 3.0) / ( pow( 0.1 * r_sunearth, 2.0) * (1.0e6 * hbar)) / (3.1709791983765E-8 * 1.0e-4) # for units of 1/(keV y m²)
  
  for e in energies:

    var iEindexx = ((e - 1.0) / 9.0).toInt 
    var diff_flux = 0.0
    
    var r_last = 0.0
    var summm = 0.0
    var sum = 0.0
    for r in 0..<df["Rho"].len:
      var n_e_keV = pow(10.0, (n_es[r].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion

      var t_keV = pow(10.0, (temperatures[r].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
      var e_keV = e * 0.001
      let r_mm = (r.float * 0.0005 + 0.0015) * r_sun
      let r_perc = (r.float * 0.0005 + 0.0015)
      if e_keV > 0.4:
        let k = sqrt((e_keV * e_keV) - ((4.0 * PI * alpha * n_e_keV) / m_e_keV)) #However, at energies near and below a typical solar plasma frequency, i.e., for energies near or below 0.3 keV,this calculation is not appropriate because the charged particles were treated as staticsources of electric fields, neglecting both recoil effects and collective motions. 
        diff_flux +=  emratesS[r][iEindexx]  * (r_perc - r_last)  *  r_perc * r_perc * e_keV * e_keV * 0.5 / (PI * PI) #k instead of e

      else : 
        diff_flux +=  emratesS[r][iEindexx]  * (r_perc - r_last)  *  r_perc * r_perc * e_keV * e_keV * 0.5 / (PI * PI)
      summm = summm + (r_perc - r_last)
      sum += (r_perc - r_last)

      r_last = r_perc 
    diff_flux = diff_flux * factor
    diff_fluxs.add(diff_flux)

    result = diff_fluxs

proc hash(x: ElementKind): Hash = 
  var h: Hash = 0
  result = h !& int(x)
  result = !$result

proc main*(energies : seq[float]) : seq[seq[float]] = 

  ## First lets access the solar model and calculate some necessary values
  const solarModel = "./ReadSolarModel/resources/AGSS09_solar_model_stripped.dat"
  var df = readSolarModel(solarModel)
  #df = df.filter(f{"Radius" <= 0.2})
  echo df.pretty(precision = 10)

  ## now let's plot radius against temperature colored by density
  ggplot(df, aes("Radius", "Temp", color = "Rho")) +
    geom_line() +
    ggtitle("Radius versus temperature of solar mode, colored by density") +
    ggsave("radius_temp_density.pdf")

  var 
    n_Z = newSeqWith(df["Rho"].len, newSeq[float](29)) #29 elements
    n_e : float
    n_e_old : float
    n_es : seq[int]
    n_eInt : int
    distNe : float
    distTemp : float
    temperature : int
    temperatures : seq[int]
    rs3 : seq[float]

  let noElement = @[3, 4, 5, 9, 15, 17, 19, 21, 22, 23, 27]
  const
    alpha = 1.0 / 137.0
    g_ae = 1e-13 # Redondo 2013: 0.511e-10 
    m_e_keV = 510.998 #keV
    e_charge = sqrt(4.0 * PI * alpha)#1.0
    kB = 1.380649e-23
    r_sun = 6.957e11 #mm
    r_sunearth = 1.5e14 #mm
    hbar = 6.582119514e-25 # in GeV * s
    keV2cm = 1.97327e-8 # cm per keV^-1
    amu = 1.6605e-24 #grams
  let atomicMass = [1.0078,4.0026,3.0160,12.0000,13.0033,14.0030,15.0001,15.9949,16.9991,17.9991,20.1797,22.9897,24.3055,26.9815,28.085,30.9737,32.0675,35.4515,39.8775,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.9380,55.845,58.9331,58.6934] #all the 29 elements from the solar model file
  let elements = ["H1", "He4","He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"]
  let charges = [1.0, 2.0, 2.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 8.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0] 

  for iRadius in 0..< df["Rho"].len:
    n_Z[iRadius][1] = (df[elements[0]][iRadius].toFloat / atomicMass[0]) * (df["Rho"][iRadius].toFloat / amu) # Hydrogen
    for iZmult in 1..3:
      if iZmult == 1:
        n_Z[iRadius][iZmult * 2] = ((df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat) / ((atomicMass[iZmult * 2 - 1] * df[elements[iZmult * 2 - 1]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)) + (atomicMass[iZmult * 2] * df[elements[iZmult * 2]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
      else: n_Z[iRadius][iZmult + 4] = ((df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat) / ((atomicMass[iZmult * 2 - 1] * df[elements[iZmult * 2 - 1]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)) + (atomicMass[iZmult * 2] * df[elements[iZmult * 2]][iRadius].toFloat / (df[elements[iZmult * 2 - 1]][iRadius].toFloat + df[elements[iZmult * 2]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu) 
    n_Z[iRadius][8] = ((df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat) / ((atomicMass[7] * df[elements[7]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[8] * df[elements[8]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)) + (atomicMass[9] * df[elements[9]][iRadius].toFloat / (df[elements[7]][iRadius].toFloat + df[elements[8]][iRadius].toFloat + df[elements[9]][iRadius].toFloat)))) * (df["Rho"][iRadius].toFloat / amu)
    for iZ in 10..<29:
      n_Z[iRadius][iZ] = (df[elements[iZ]][iRadius].toFloat / atomicMass[iZ]) * (df["Rho"][iRadius].toFloat / amu) # The rest
    n_e = 0.0
    for Z in 0..<elements.len:
      n_e += (df["Rho"][iRadius].toFloat/amu) * charges[Z] * df[elements[Z]][iRadius].toFloat / atomicMass[Z] # (g/cm³ /g) = 1/cm³
    n_e_old = (df["Rho"][iRadius].toFloat/amu) * (1 + df[elements[0]][iRadius].toFloat/2)
    for iTemp in 0..90:
      distTemp = log(df["Temp"][iRadius].toFloat, 10.0) / 0.025 - float(140 + 2 * iTemp)
      if abs(distTemp) <= 1.0: 
        temperature = 140 + 2 * iTemp
    temperatures.add(temperature)
    for iNe in 0..17:
      distNe = log(n_e, 10.0) / 0.25 - float(74 + iNe * 2) 
      if abs(distNe) <= 1.0: 
        n_eInt = 74 + iNe * 2 
    n_es.add(n_eInt)
    rs3.add(iRadius.float * 0.0005 + 0.0015)

  var dfTemp = seqsToDf({ "Radius" : rs3,
                        "Temp" : temperatures,
                        "Ne" : n_es})
  echo dfTemp

  ggplot(dfTemp, aes("Radius", "Ne", color = "Temp")) +
    geom_point() +
    ggtitle("Radius versus temperature of solar mode, colored by density") +
    ggsave("radius_temp_ne.pdf")



  var densities: HashSet[int]
  var opElements: Table[ElementKind, Table[int, OpacityFile]]
  var opElNew: Table[ZTempDensity, OpacityFile]


  for temp in toSet(temperatures):
    for (Z_str, Z) in iterEnum(ElementKind):
      let testF = &"./OPCD_3.3/mono/fm{Z:02}.{temp}"
      if existsFile(testF):
        let opFile = parseOpacityFile(testF, kind = ofkOriginal)
        for k in keys(opFile.densityTab):
          densities.incl k
        let zKind = ElementKind(Z)
        if zKind notin opElements:
          opElements[zKind] = initTable[int, OpacityFile]()
        opElements[zKind][temp] = opFile
      for ne in toSet(n_es):
        let opFile = &"./OPCD_3.3/OP/opacity_table_{Z_str[1 .. ^1]}_{temp}_{ne}.dat"
        if existsFile(opFile):
          opElNew[(Z, temp, ne)] = parseOpacityFile(opFile, kind = ofkNew)
          
        


  ## Calculate the absorbtion coefficients depending on the energy and the radius out of the opacity values


  var absCoefs = newSeqWith(df["Rho"].len, newSeq[float](1112)) #29 elements
  var emratesS = newSeqWith(df["Rho"].len, newSeq[float](1112))
  var term1s = newSeqWith(df["Rho"].len, newSeq[float](1112))
  var comptons = newSeqWith(df["Rho"].len, newSeq[float](1112))
  var term3s = newSeqWith(df["Rho"].len, newSeq[float](1112))
  var ffterms = newSeqWith(df["Rho"].len, newSeq[float](1112))

  var posOP = newSeqWith(df["Rho"].len, newSeq[int](1112))
  var ironOpE : seq[float]
  var n_e_keV : float
  echo df["Rho"].len
  var zs : seq[int]
  var rs : seq[float]
  var rs2 : seq[float]
  var nZs : seq[float]
  var engs : seq[float]
  var ops : seq[float]

  for R in 0..<df["Rho"].len:
    n_eInt = n_es[R]
    temperature = temperatures[R]
    for iE in energies:
      var sum = 0.0
      var absCoef = 0.0
      n_e_keV = pow(10.0, (n_es[R].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion
      var temp_keV = pow(10.0, (temperatures[R].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
      var temp_K = pow(10.0, (temperatures[R].toFloat * 0.025))
      
      var energy_keV = iE * 0.001 #w * temp_keV
      var energy_J = 1.60218e-16 * energy_keV
      var w = energy_keV / temp_keV #toFloat(dfMesh["u"][iE.int])
      var iEindex = ((iE - 1.0) / 9.0).toInt 
      var table = 1.0
      var prevMesh = 1.0      
      let debye_scale = sqrt( (4.0 * PI * alpha / temp_keV) * (n_e_keV + n_Z[R][1] * 7.645e-24 + 4.0 * n_Z[R][2] * 7.645e-24 )) 
      
      let y = debye_scale / (sqrt( 2.0 * m_e_keV * temp_keV))
      if w >= 20.0 or w <= 0.0732: #because the tables dont go beyond that, apparently because the axion production beyond that is irrelevant #except for He, maybe find a better solution
        for (Z_str, Z) in iterEnum(ElementKind):
          if Z in noElement:
            continue
          sum = sum + n_Z[R][Z] * 0.0       
        absCoefs[R][iEindex] = sum * 1.97327e-8 * 0.528e-8 * 0.528e-8 * (1.0 - exp(-energy_keV / temp_keV))
      else : 

        for (Z_str, Z) in iterEnum(ElementKind):
          if Z in noElement:
            continue
          table = spline.eval(w)
          if iEindex == 11 :
            zs.add(Z)
            rs.add((R.float * 0.0005 + 0.0015))
            nZs.add(n_Z[R][Z])
          

          if Z == 1 and (iEindex == 12 or iEindex == 68 or iEindex == 140 or iEindex == 239 or iEindex == 599): #and (iEindex == 12 or iEindex == 68) :
            rs2.add((R.float * 0.0005 + 0.0015))
            engs.add(energy_keV)
            ops.add(opElNew[(Z, temperature, n_eInt)].densityOp.interp.eval(energy_keV))
          let opacityL = opElNew[(Z, temperature, n_eInt)].densityOp.interp.eval(energy_keV)
          let opacity = opElements[ElementKind(Z)][temperature].densityTab[n_eInt].interp.eval(table) 


          #var opacity_cm = opacity  # correct conversion
          # opacities in atomic unit for lenth squared: 0.528 x10-8cm * 0.528 x10-8cm = a0² # 1 m = 1/1.239841336215e-9 1/keV and a0 = 0.528 x10-10m
          if Z > 2:
            sum +=  n_Z[R][Z] * opacity
          

        
        absCoef = sum * 1.97327e-8 * 0.528e-8 * 0.528e-8 * (1.0 - exp(-energy_keV / temp_keV)) # is in keV  
        
        absCoefs[R][iEindex] = absCoef
      
      ## Now it's left to calculate the emission rates
      ## making the same approximation as for n_e calculation 

      ##ion density weighted by charge^2 from Raffelt
      let nZZ2_raffelt_keV = (df["Rho"][R].toFloat / amu) * 7.683e-24 #seems to be correct
      let ffterm = freefreeEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV, nZZ2_raffelt_keV, w, y) 
      ## includes contribution from ff, fb and bb processes and a part of the Comption contribution ## keV³ / keV² = keV :
      let term1 = term1(g_ae, energy_keV, (absCoefs[R][iEindex]), e_charge, m_e_keV, temp_keV)  
      let term2 = term2(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV)# completes the Compton contribution #keV
      let term3 = bremsEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV, w, y) # contribution from ee-bremsstahlung
      let compton = comptonEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV) 
      term1s[R][iEindex] = term1
      comptons[R][iEindex] = compton
      term3s[R][iEindex] = term3
      ffterms[R][iEindex] = ffterm 
      let total_emrate = compton +  term1 + term3 + ffterm#) keV 
      let total_emrate_s = total_emrate #/ (6.58e-19) # in 1/sec 
      emratesS[R][iEindex] = total_emrate_s
      # if want to have absorbtion coefficient of a radius and energy: R = (r (in % of sunR) - 0.0015) / 0.0005
      
      #if iEindex == 110:
        #echo fNew(w, y), " for r ", (R.float * 0.0005 + 0.0015), " and e ", energy_keV 
      #echo total_emrate


  
  var dfNZ = seqsToDf({ "Radius" : rs,
                        "nZ" : nZs,
                        "Z" : zs})
  echo dfNZ
  ggplot(dfNZ, aes("Radius", "nZ", color = "Z")) +
    geom_point() +
    ggtitle("Radius versus atomic density for different Z") +
    ggsave("radius_nZ_Z.pdf") 

  var dfOp = seqsToDf({ "Radius" : rs2,
                        "opacity" : ops,
                        "energies" : engs})
  echo dfOp
  ggplot(dfOp, aes("Radius", "opacity", color = "energies")) +
    geom_point() +
    ggtitle("Radius versus opacity for different energies") +
    ggsave("radius_op_energies.pdf")   

  let diff_fluxs = getFluxFraction(energies, df, n_es, temperatures, emratesS)
  let term1flux = getFluxFraction(energies, df, n_es, temperatures, term1s)
  let comptonflux = getFluxFraction(energies, df, n_es, temperatures, comptons)
  let term3flux = getFluxFraction(energies, df, n_es, temperatures, term3s)
  let fftermflux = getFluxFraction(energies, df, n_es, temperatures, ffterms)
  var energieslong : seq[float]
  var fluxes : seq[float]
  var kinds : seq[string]
  for E in energies:
    var e_keV = E * 0.001
    var iEindex = ((E - 1.0) / 9.0).toInt
    energieslong.add(e_keV)
    fluxes.add(diff_fluxs[iEindex])
    kinds.add("Total Flux")
    energieslong.add(e_keV)
    fluxes.add(term1flux[iEindex])
    kinds.add("FB BB Flux")
    energieslong.add(e_keV)
    fluxes.add(comptonflux[iEindex])
    kinds.add("Compton Flux")
    energieslong.add(e_keV)
    fluxes.add(term3flux[iEindex])
    kinds.add("EE Flux")
    energieslong.add(e_keV)
    fluxes.add(fftermflux[iEindex])
    kinds.add("FF Flux")



  let dfEmrate = seqsToDf({ "energy" : energies,
                            "emrate" : emratesS[200] })
  ggplot(dfEmrate, aes("energy", "emrate")) +
    geom_line() +
    ggsave("emrate_R10.pdf")

  let dfAbscoef = seqsToDf({ "energy" : energies,
                            "absCoefs" : absCoefs[10] })
  ggplot(dfAbscoef, aes("energy", "absCoefs")) +
    geom_line() +
    ggsave("abscoefs_R10.pdf")


  let dfDiffflux = seqsToDf({ "Axion energy [eV]" : energieslong,
                              "Fluxfraction [keV⁻¹y⁻¹m⁻²]" : fluxes,
                              "type" : kinds })
  ggplot(dfDiffflux, aes("Axion energy [eV]", "Fluxfraction [keV⁻¹y⁻¹m⁻²]", color = "type")) +
    geom_line() +
    ggtitle("The flux fraction of axions from the sun") +
    ggsave("diffFlux.pdf")
  
  result = emratesS

#echo main()