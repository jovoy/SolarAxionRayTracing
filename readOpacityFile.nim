import std / [strutils, math, tables, sequtils, strformat, hashes, macros, os, strscans]
import pkg / [polynumeric, ggplotnim, numericalnim, glm]

import seqmath except linspace
import arraymancer except readCsv, linspace
import json except `{}` # to not get into trouble with datamancer `f{}`

# for OPCD data file parsing
import cligen/[osUt, mfile, mslice]

# local import
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
    interp: InterpolatorType[float]

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

const atomicMass =
  [1.0078, 4.0026, 3.0160, 12.0000, 13.0033, 14.0030, 15.0001, 15.9949, 16.9991,
   17.9991, 20.1797, 22.9897, 24.3055, 26.9815, 28.085, 30.9737, 32.0675,
   35.4515, 39.8775, 39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961,
   54.9380, 55.845, 58.9331, 58.6934] #all the 29 elements from the solar model file
const elements =
  ["H1", "He4","He3", "C12", "C13", "N14", "N15", "O16", "O17", "O18", "Ne",
   "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V",
   "Cr", "Mn", "Fe", "Co", "Ni"]
const charges =
  [1.0, 2.0, 2.0, 6.0, 6.0, 7.0, 7.0, 8.0, 8.0, 8.0, 10.0, 11.0, 12.0, 13.0,
   14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
   26.0, 27.0, 28.0]

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

proc parseTableLine(energy, opacity: var float, ms: MSlice) {.inline.} =
  ## parses the energy and opacity float values from `line` into `energy`
  ## (1st col) and `opacity` (2nd col) using `cligen.parseFloat` and `mslices`
  var
    i = 0 # manual counter for processed numbers
    x1, x2: int
  for el in mslices(ms, sep = ' '):
    if el.len == 0: continue # ignore empty slice
    if i == 0: opacity = parseFloat(el)
    elif i == 1:
      energy = parseFloat(el)
      swap(energy, opacity) # swap both because if energy present it's the first column
    else:
      raise newException(ValueError, "Parsing opacity table in line " & $ms & " failed!")
    inc i

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

proc parseDensityTab(mf: var MSlice,
                     temp: int,
                     kind: OpacityFileKind,
                     tableCount = 1000): DensityOpacity =
  # now parse the table according to table count
  result = DensityOpacity(energies: newSeqOfCap[float](tableCount),
                          opacities: newSeqOfCap[float](tableCount))
  var
    buf = newString(200)
    energy: float
    opac: float
    lineCnt = 0
    ms: MSlice
  while mf.len > 0:
    discard mf.nextSlice(ms)
    parseTableLine(energy, opac, ms)
    case kind
    of ofkOriginal:
      if tableCount == 10000:
        # set energy manually to `j`, since we simply have 1 eV steps
        energy = float lineCnt + 1
      # else the input file contains the energy
    of ofkNew:
      energy = convertLogUToE(energy, temp.float)
    result.energies.add energy
    result.opacities.add opac
    inc lineCnt
    if kind == ofkOriginal and lineCnt == tableCount:
      break
  # finalize densityOpacity by creating spline and adding to result
  result.interp = newLinear1D(result.energies,
                              result.opacities)

proc parseOpacityNew(mfile: MFile,
                     fname: string,
                     kind: static OpacityFileKind): OpacityFile =
  var buf = newString(200)
  doAssert false, "Currently unsupported."
  when false:
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

proc parseOpacityOriginal(mfile: MFile,
                          fname: string,
                          kind: static OpacityFileKind): OpacityFile =
  let temp = parseInt(fname[5 .. ^1])
  result = OpacityFile(fname: fname,
                       kind: kind,
                       element: ElementKind(fname[2 .. 3].parseInt),
                       temp: temp) #pow(10.0, parseFloat(fname[5 .. ^1]) / 40.0))
  echo "Parsing file : ", fname
  var ms: MSlice
  var ds = mfile.toMSlice()
  var first = true
  while ds.len > 0:
    if first:
      # skip file header, drop the slice
      discard ds.nextSlice(ms)
      first = false
    # read table header, 3 lines
    discard ds.nextSlice(ms)
    let h1 = parseTableHeader($ms, H1) # copy to nim string irrelevant for perf
    discard ds.nextSlice(ms)         # due to multiple orders of magnitude more
    let h2 = parseTableHeader($ms, H2) # data than header
    discard ds.nextSlice(ms)
    var tableCount = ($ms).strip.parseInt # inefficient, but rare
    tableCount = if tableCount == 0: 10000 else: tableCount
    result.densityTab[h1.density] = ds.parseDensityTab(temp, kind, tableCount)

proc parseOpacityFile(path: string, kind: OpacityFileKind): OpacityFile =
  ## we parse the monochromatic opacity file using strscans
  ## - first we drop the first line as the file header. Information in this?
  ## - then read table header (3 lines)
  ## - then num lines
  let ds = mopen(path)
  if ds.mem.isNil:
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
  var df = readCsv(fname, sep = ' ')
  doAssert df.len == 10001
  result = df["u"].toTensor(float).toRawSeq


let lineNumbers = linspace(0.0, 10000.0, 10001)
const meshFile = "./OPCD_3.3/mono/fm01.mesh"
## NOTE: all the `fm??.mesh` files are identical!
## check `tools/diff_files.nim` for proof.
let dfMesh = readMeshFile(meshFile)

let spline = newLinear1D(dfMesh, lineNumbers)

template inner_integral(t: float, y: float): float =
  (1.0/2.0) * ( ((y * y) / (t * t + y * y)) + ln( t * t + y * y ) )

proc outer(x, w, y: float): float {.inline.} =
  let coeff = x * exp(-x * x)
  # wrap the inner call to have an `IntegrateFunc` kind
  #let fn = proc(x: float, optional: seq[float]): float =
    #result = inner(x, y)
  let
    frm = sqrt(x * x + w) - x
    to = sqrt(x * x + w) + x
  let integral = inner_integral(to, y) - inner_integral(frm, y)#adaptiveGauss(fn, frm, to) #check if analytical integral gives the same results: difference of e-15

  result = coeff * integral

proc fNew(w: float, y: float): float =
  # integrate `outer` frm `0` to `inf`
  # rewrite via
  # int_a^infty dx f(x) = int_0^1 f(a + (1 - t) / t) / t^2
  # in our case a = 0
  # so express by wrapping `outer` in a new proc
  let fnToInt = proc(t: float, optional: NumContext[float, float]): float =
    #echo t
    if t != 0:
      result = outer(x = (1 - t) / t, w = w, y = y) / (t * t)
    else:
      # workaround singularity by adding epsilon
      result = outer(x = (1 - t) / (t + 1e-8), w = w, y = y) / (t * t)
  # and integrate that from 0 to 1
  result = adaptiveGauss(fnToInt, 0.0, 1.0)

proc bfield(r:float): float = #r in sun radius percentage
  let
    radius_cz = 0.712 # start of convective zone = end of radiative zone
    size_tach = 0.02 # size of convective/radiative transistion region
    radius_outer = 0.96 # upper layers of the Sun
    size_outer = 0.035 # size of upper layers of the Sun
    # B-field reference values
    bfield_rad_T = 3.0e3
    bfield_tach_T = 50.0
    bfield_outer_T = 4.0
    lambda1 = 10.0*radius_cz + 1.0
    lambda_factor = (1.0 + lambda1)*pow(1.0 + 1.0/lambda1, lambda1);
  var bfield = 0.0
  if r < (radius_cz + size_tach):
    var x = pow(r/radius_cz, 2.0)
    if x < 1.0: bfield = bfield_rad_T * lambda_factor * x * pow(1.0-x, lambda1)
    var y = pow(((r - radius_cz) / size_tach), 2.0)
    if y < 1.0: bfield = bfield_tach_T * (1.0-y)
  else:
    var z = pow((r - radius_outer)/size_outer, 2.0)
    if z < 1.0: bfield = bfield_outer_T*(1.0 - z)
    else: bfield = 0.0

  result = bfield/(1.0e6 * 1.4440271 * 1.0e-3 * sqrt(4.0 * PI)) # in keV^2



proc omegaPlasmonSq(alpha:float, ne: float, me: float): float = #Routine to return the plasma freqeuency squared (in keV^2) of the zone around the distance r from the centre of the Sun.
  result = 4.0 * alpha * PI * ne / me

## The functions for all the parts of the emission rate ##

proc comptonEmrate(alpha, gae, energy, ne, me, temp: float): float =
  result = (alpha * gae * gae * energy * energy * ne) /
           (3.0 * pow(me, 4) * (exp(energy / temp) - 1.0))

proc bremsEmrate(alpha, gae, energy, ne, me, temp, w, y: float): float =
  result = (alpha * alpha * gae * gae * 4.0 * sqrt(PI) * ne * ne *
              exp(- energy / temp) * fNew(w, sqrt(2.0) * y)) /
           (3.0 * sqrt(temp) * pow(me, 3.5) * energy)

proc term1(gae, energy, abscoef, echarge, me, temp: float): float =
  result = (gae * gae * energy * energy * abscoef) /
           (2.0 * echarge * echarge * me * me * (exp(energy / temp) - 1.0))

proc term2(alpha, gae, energy, ne, me, temp: float): float =
  result = ((exp(energy / temp) - 2.0) *
               comptonEmrate(alpha, gae, energy, ne, me, temp)) /
           (2.0 * (exp(energy / temp) - 1.0))

proc freefreeEmrate(alpha, gae, energy, ne, me, temp, nzZ2, w, y: float): float =
  result = (fNew(w, y) * alpha * alpha * gae * gae  * 8.0 * sqrt(PI) * ne * nzZ2 *
                exp(-energy/temp)) /
           (3.0 * sqrt(2.0 * temp) * pow(me, 3.5) * energy)

# Some auxilliary functions for Primakoff rate and degeneracy calculation
proc primakoff_bracket(t: float, u: float): float =
  var analytical_integral = 0.0
  if u > 1.0:  analytical_integral += (u * u - 1.0)* ln((u - 1.0)/( u + 1.0))
  var v = u + t
  if v > 1.0: analytical_integral -= (v * v - 1.0) * ln((v - 1.0) / (v + 1.0))
  analytical_integral *= 0.5 / t
  analytical_integral -= 1.0
  result = analytical_integral


proc primakoff(temp, energy, gagamma, ks2, alpha, ne, me: float, n_Z2: float, n_Z1: float): float =
  let
    prefactor6 = gagamma * gagamma * 1e-12 * alpha  / 8.0
    omPlSq = omegaPlasmonSq(alpha, ne, me)
    z = energy / temp
    om2 = energy * energy
    x = om2/omPlSq

  ## from Raffelt 2006 dont know if this is the newest
  if x < 1.0 or energy == 0.0:
    result = 0.0
  else:
    let
      phase_factor = 2.0 / (sqrt(1.0 - 1.0 / x) * (exp(z) - 1.0))
      n_dens = ne + n_Z1 * 7.645e-24 + 4.0 * n_Z2 * 7.645e-24# ne + n_Z2 * pow(197.327053e-10,3.0)  #;avg_degeneracy_factor(r)*
      s = 2.0*energy*sqrt(om2 - omPlSq)
      t = ks2 / s
      u = (2.0 * om2 - omPlSq) / s
      analytical_integral = primakoff_bracket(t, u)
    result = prefactor6 * phase_factor * n_dens * analytical_integral
  #[result = (gagamma * gagamma * 1e-12 * temp * ks2) /
          (32.0 * PI) *
          ((1.0 + (ks2 / (4.0 * energy * energy))) *
              ln(1.0 + (4.0 * energy * energy) / ks2) - 1.0) *
          2.0 * sqrt(1.0 - (4.0 * PI * ne) / (me * energy * energy)) /
          (exp(energy/temp) - 1.0)]##[]#

proc longPlasmon(energy: float, ne: float, me: float, alpha: float, bfieldR: float, temp: float, opacity: float, gagamma: float): float = #https://arxiv.org/pdf/2101.08789.pdf
  let
    omPlSq = omegaPlasmonSq(alpha, ne, me)
    prefactor = gagamma * gagamma * 1e-12
    om2 = energy * energy
    z = energy/temp
  var gammaL = (1.0 - exp(-z)) * opacity
  gammaL = max(gammaL, 1e-4) # to avoid numerical issues from very narrow resonances
  let
    xi2 = gammaL*energy
    fwhm = sqrt(om2 + xi2) - sqrt(om2 - xi2) # FWHM of Lorentz/Cauchy peak
  #### if (gsl_pow_2(om2 - om_pl_sq) > 100.0 * om2*gammaL*gammaL) { return 0; } //just integrate around resonance
  if abs(energy - sqrt(omPlSq)) > 18.0*fwhm: return 0 # Just integrate around resonance
  let
    average_bfield_sq = bfieldR*bfieldR/3.0
    fraction = energy*xi2 / ( pow(om2 - omPlSq, 2.0) + xi2*xi2 )
  result = prefactor * average_bfield_sq * fraction / (exp(z) - 1.0)

proc transPlasmon(energy: float, ne: float, me: float, alpha: float, bfieldR: float, temp: float, opacity: float, gagamma: float): float =
  let
    geom_factor = 1.0 # factor accounting for observers position (1.0 = angular average)
    photon_polarization = 2.0
    omPlSq = omegaPlasmonSq(alpha, ne, me)
  if omPlSq > energy*energy: return 0 # energy can't be lower than plasma frequency
  let
    u = energy/temp
    gamma = (1.0 - exp(-u))*opacity
    deltaPsq = energy*energy * pow(sqrt(1.0-omPlSq/(energy*energy))-1.0, 2.0)  # transfered momentum squared
    #deltaPsq = pow(0.5*omPlSq/energy, 2.0) #transfered momentum squared
    average_b_field_sq = pow(bfieldR, 2.0) / 3.0
    deltaTsq = gagamma*gagamma * 1e-12 * average_b_field_sq / 4.0
  result = geom_factor * photon_polarization * gamma * deltaTsq / ( (deltaPsq+pow(0.5*gamma, 2.0)) * (exp(u) - 1.0) )

proc iron(ganuclei: float, temp: float, energy: float, rho: float): float = #https://iopscience.iop.org/article/10.1088/1475-7516/2009/12/002/pdf
  let
    tau_gamma = 1.3e-6 * 1.519e18 #1/keV was s
    n = 3.0e17 * 1.7826e-30 #1/keV was g⁻1
    e_gamma = 14.4 #keV
    m_Fe = 56.9353928 * 1.6605e-24 * 5.60958616722e29 #was g now keV
    u = e_gamma/temp #what about k? is correct no k
    w_1 = 4.0 * exp(-u)/(2.0 + 4.0 * exp(-u)) # the former is used in the paper as an approximation #2.0 * exp(-u) #
    gamma_frac = 1.82 * ganuclei * ganuclei
    sigma = e_gamma * sqrt(temp / m_Fe)  #in keV
    n_a = n * w_1 * gamma_frac / tau_gamma  #in no units #/ 5.60958616722e29 * 4.135665538536e-18
  #echo w_1, " ", 4.0 * exp(-u)/(2.0 + 4.0 * exp(-u))
  result = n_a * exp(- pow(energy - e_gamma, 2.0) / (2.0 * sigma * sigma)) * rho * sqrt(2.0 * PI) * PI/ (sigma * energy * energy) #because of the way the flux is portrait


proc getFluxFraction(energies: seq[float], df: DataFrame,
                     n_es, temperatures: seq[int],
                     emratesS: Tensor[float],
                     typ: string = ""): DataFrame =
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
  let factor = pow(r_sun * 0.1 / (keV2cm), 3.0) /
               (pow(0.1 * r_sunearth, 2.0) * (1.0e6 * hbar)) /
               (3.1709791983765E-8 * 1.0e-4) # for units of 1/(keV y m²)
  var diff_fluxs: seq[float]
  var radii: seq[float]
  var E: seq[float]
  for e in energies:
    var iEindexx = ((e - 1.0)).toInt
    var diff_flux = 0.0
    var r_last = 0.0
    var summm = 0.0
    var sum = 0.0
    for r in 0 ..< df["Rho"].len:
      let
        n_e_keV = pow(10.0, (n_es[r].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion
        t_keV = pow(10.0, (temperatures[r].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
        e_keV = e * 0.001
        r_mm = (r.float * 0.0005 + 0.0015) * r_sun
        r_perc = (r.float * 0.0005 + 0.0015)
      if e_keV > 0.4:
        # However, at energies near and below a typical solar plasma frequency,
        # i.e., for energies near or below 0.3 keV,this calculation is not
        # appropriate because the charged particles were treated as static
        # sources of electric fields, neglecting both recoil effects and collective motions.
        ## TODO: what the heck is this? `k` is nowhere to be seen
        let k = sqrt((e_keV * e_keV) - ((4.0 * PI * alpha * n_e_keV) / m_e_keV))
        diff_flux = emratesS[r, iEindexx] * (r_perc - r_last) * r_perc * r_perc *
                     e_keV * e_keV * 0.5 / (PI * PI) #k instead of e
      else :
        diff_flux = emratesS[r, iEindexx] * (r_perc - r_last) * r_perc * r_perc *
                     e_keV * e_keV * 0.5 / (PI * PI)
      summm = summm + (r_perc - r_last)
      sum += (r_perc - r_last)
      r_last = r_perc

      diff_fluxs.add diffFlux * factor
      radii.add r.float
      E.add e
    #diff_flux = diff_flux * factor
    #diff_fluxs.add(diff_flux)

    #result = diff_fluxs
  result = seqsToDf({"diffFlux" : diffFluxs, "Radius" : radii, "Energy" : E })
  result["type"] = constantColumn(typ, result.len)

proc getFluxFractionR(energies: seq[float], df: DataFrame,
                     n_es, temperatures: seq[int],
                     emratesS: Tensor[float],
                     typ: string = ""): DataFrame =
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
  let factor = pow(r_sun * 0.1 / (keV2cm), 3.0) /
               (pow(0.1 * r_sunearth, 2.0) * (1.0e6 * hbar)) /
               (3.1709791983765E-8 * 1.0e-4) # for units of 1/(keV y m²)
  var diff_fluxs: seq[float]
  var radii: seq[float]
  var E: seq[float]
  for (idx, e) in pairs(energies):
    var
      diff_flux = 0.0
      diff_fluxR = 0.0
      r_last = 0.0
      summm = 0.0
      sum = 0.0
    for r in 0 ..< df["Rho"].len:
      let
        n_e_keV = pow(10.0, (n_es[r].toFloat * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion
        t_keV = pow(10.0, (temperatures[r].toFloat * 0.025)) * 8.617e-8 # was K # correct conversion
        e_keV = e * 0.001
        r_mm = (r.float * 0.0005 + 0.0015) * r_sun
        r_perc = (r.float * 0.0005 + 0.0015)
      if e_keV > 0.4:
        # However, at energies near and below a typical solar plasma frequency,
        # i.e., for energies near or below 0.3 keV,this calculation is not
        # appropriate because the charged particles were treated as static
        # sources of electric fields, neglecting both recoil effects and collective motions.
        ## TODO: what the heck is this? `k` is nowhere to be seen
        let k = sqrt((e_keV * e_keV) - ((4.0 * PI * alpha * n_e_keV) / m_e_keV))
        diff_flux = emratesS[r, idx] * (r_perc - r_last) * r_perc * r_perc *
                     e_keV * e_keV * 0.5 / (PI * PI) #k instead of e
      else :
        diff_flux = emratesS[r, idx] * (r_perc - r_last) * r_perc * r_perc *
                     e_keV * e_keV * 0.5 / (PI * PI)
      summm = summm + (r_perc - r_last)
      sum += (r_perc - r_last)
      diffFluxR += diff_flux #* (r_perc - r_last) * r_sun
      r_last = r_perc

    diff_fluxs.add diffFluxR * factor

    E.add e
    #diff_flux = diff_flux * factor
    #diff_fluxs.add(diff_flux)

    #result = diff_fluxs
  result = seqsToDf({"diffFlux" : diffFluxs, "Energy" : E })
  result["type"] = typ

#proc getFluxFraction(energies: seq[float], df: DataFrame,
#                     n_es, temperatures: seq[int],
#                     emratesS: Tensor[float]): seq[float] =
#  result = df.group_by("Radius").summarize(f{sum(`diffFlux`)})["diffFlux", float].toRawSeq

proc hash(x: ElementKind): Hash =
  var h: Hash = 0
  result = h !& int(x)
  result = !$result

proc main*(): Tensor[float] =

  ## First lets access the solar model and calculate some necessary values
  const solarModel = "./ReadSolarModel/resources/AGSS09_solar_model_stripped.dat"
  var df = readSolarModel(solarModel)
  let nElems = 1500 # TODO: clarify exact number
  let energies = linspace(1e-3, 15.0, nElems)

  let nRadius = df["Rho"].len

  ## now let's plot radius against temperature colored by density
  ggplot(df, aes("Radius", "Temp", color = "Rho")) +
    geom_line() +
    ggtitle("Radius versus temperature of solar mode, colored by density") +
    ggsave("out/radius_temp_density.pdf")

  var
    n_Z = newSeqWith(nRadius, newSeq[float](29)) #29 elements
    n_e: float
    n_e_old: float
    n_esfloat: seq[float]
    n_es: seq[int]
    n_eInt: int
    distNe: float
    distTemp: float
    temperature: int
    temperatures: seq[int]
    tempFloat: seq[float]
    rs3: seq[float]
    alphaR: seq[float]
    bfields: seq[float]
    ompls: seq[float]

  let noElement = @[3, 4, 5, 9, 15, 17, 19, 21, 22, 23, 27]
  const
    alpha = 1.0 / 137.0
    g_ae = 1e-13 # Redondo 2013: 0.511e-10  #1e-11 #
    gagamma = 1e-12 #the latter for DFSZ  #1e-9 #5e-10 #
    m_a = 0.0853 #eV
    ganuclei = 1e-15 #1.475e-8 * m_a #KSVZ model #no units  #1e-7
    m_e_keV = 510.998 #keV
    e_charge = sqrt(4.0 * PI * alpha)#1.0
    kB = 1.380649e-23
    r_sun = 6.957e11 #mm
    r_sunearth = 1.5e14 #mm
    hbar = 6.582119514e-25 # in GeV * s
    keV2cm = 1.97327e-8 # cm per keV^-1
    amu = 1.6605e-24 #grams
  # send halp

  echo "Walking all radii"
  let rho = df["Rho"].toTensor(float)

  for iRadius in 0..< rho.size:
    template eAt(arg: untyped): untyped = df[elements[arg]][iRadius, float]
    template aAt(arg: untyped): untyped = atomicMass[arg]
    template n(idx: int): untyped =
      (eAt(idx) / aAt(idx)) * (rho[iRadius] / amu)

    n_Z[iRadius][1] = n(0) # Hydrogen
    for iZmult in 1..3:
      let iz = iZmult * 2
      let nVal = (eAt(iz - 1) + eAt(iz)) /
                 ((aAt(iz - 1) * eAt(iz - 1) + aAt(iz) * eAt(iz)) /
                  (eAt(iz - 1) + eAt(iz))) *
                 rho[iRadius] / amu
      if iZmult == 1:
        n_Z[iRadius][iz] = nVal
      else:
        n_Z[iRadius][iZmult + 4] = nVal

    n_Z[iRadius][8] = (eAt(7) + eAt(8) + eAt(9)) /
                       ((eAt(7) * aAt(7) + eAt(8) * aAt(8) + eAt(9) * aAt(9)) /
                        (eAt(7) + eAt(8) + eAt(9))) *
                      rho[iRadius] / amu
    for iZ in 10..<29:
      n_Z[iRadius][iZ] = n(iZ)
    n_e = 0.0
    for Z in 0..<elements.len:
      n_e += (rho[iRadius]/amu) * charges[Z] * df[elements[Z]][iRadius, float] / atomicMass[Z] # (g/cm³ /g) = 1/cm³
    n_e_old = (rho[iRadius]/amu) * (1 + df[elements[0]][iRadius, float]/2)
    n_esfloat.add(n_e)
    for iTemp in 0..90:
      distTemp = log(df["Temp"][iRadius, float], 10.0) / 0.025 - float(140 + 2 * iTemp)
      if abs(distTemp) <= 1.0:
        temperature = 140 + 2 * iTemp
    temperatures.add(temperature)
    tempFloat.add(df["Temp"][iRadius, float])
    for iNe in 0..17:
      distNe = log(n_e, 10.0) / 0.25 - float(74 + iNe * 2)
      if abs(distNe) <= 1.0:
        n_eInt = 74 + iNe * 2
    n_es.add(n_eInt)
    rs3.add(iRadius.float * 0.0005 + 0.0015)
    ompls.add(omegaPlasmonSq(alpha, (n_e * 7.683e-24), m_e_keV))
    bfields.add(bfield(0.0015 + iRadius.float * 0.0005))
    var metallicity = 1.0 - n_Z[iRadius][1] - n_Z[iRadius][2]
    var alphaRs = 0.0
    #echo n_e * 7.683e-24, " ", pow(10.0, (n_eInt.float * 0.25)) * 7.683e-24
    #for k in 2..< 29:
      #alphaRs += n_Z[iRadius][k+1] / metallicity #* ionisationsqr_element(r, element) / aAt(k)
    #alphaR.add(alphaRs)

  #[var dfTemp = seqsToDf({ "Radius": rs3,
                        "Temp": temperatures,
                        "Ne": n_es})
  echo dfTemp

  ggplot(dfTemp, aes("Radius", "Ne", color = "Temp")) +
    geom_point() +
    ggtitle("Radius versus temperature of solar mode, colored by density") +
    ggsave("out/radius_temp_ne.pdf")

  var dfPlas = seqsToDf({ "Radius": rs3,
                        "Omega": ompls,
                        "B": bfields})
  ggplot(dfPlas, aes("Radius", "Omega")) +
    geom_line() +
    scale_y_log10()+
    ggtitle("Radius versus plasma frequency") +
    ggsave("out/radius_omegapl.pdf")
  ggplot(dfPlas, aes("Radius", "B")) +
    geom_line() +
    ggtitle("Radius versus B field") +
    ggsave("out/radius_B.pdf")]#



  var densities: HashSet[int]
  var opElements = newTable[int, TableRef[ElementKind, OpacityFile]]()
  echo "Walking all temps..."
  for temp in toSet(temperatures):
    if temp notin opElements:
      opElements[temp] = newTable[ElementKind, OpacityFile]() #initTable[int, OpacityFile]()
    for (Z_str, Z) in iterEnum(ElementKind):
      let testF = &"./OPCD_3.3/mono/fm{Z:02}.{temp}"
      if existsFile(testF):
        let opFile = parseOpacityFile(testF, kind = ofkOriginal)
        for k in keys(opFile.densityTab):
          densities.incl k
        let zKind = ElementKind(Z)
        opElements[temp][zKind] = opFile
      #for ne in toSet(n_es):
      #  let opFile = &"./OPCD_3.3/OP/opacity_table_{Z_str[1 .. ^1]}_{temp}_{ne}.dat"
      #  if existsFile(opFile):
      #    opElNew[(Z, temp, ne)] = parseOpacityFile(opFile, kind = ofkNew)
  ## Calculate the absorbtion coefficients depending on the energy and the radius out of the opacity values
  var
    absCoefs = zeros[float](nRadius, nElems) #29 elements
    emratesS = zeros[float](nRadius, nElems)
    emratesInS = zeros[float](nRadius, nElems)
    term1s = zeros[float](nRadius, nElems)
    comptons = zeros[float](nRadius, nElems)
    term3s = zeros[float](nRadius, nElems)
    ffterms = zeros[float](nRadius, nElems)
    primakoffs = zeros[float](nRadius, nElems)
    longPlasmons = zeros[float](nRadius, nElems)
    transPlasmons = zeros[float](nRadius, nElems)
    iron57s = zeros[float](nRadius, nElems)
    abc = zeros[float](nRadius, nElems)
    posOP = zeros[int](nRadius, nElems)
    ironOpE: seq[float]
    n_e_keV: float
    temp_keV: float
    totalRadiiFlux = newSeq[float](nRadius)
    radiiFlux: float
    radii = newSeq[float](nRadius)
    zs = newSeqOfCap[int](100_000)
    rs = newSeqOfCap[float](100_000)
    rs2 = newSeqOfCap[float](100_000)
    nZs = newSeqOfCap[float](100_000)
    engs = newSeqOfCap[float](100_000)
    ops = newSeqOfCap[float](100_000)

  echo "Walking all radii again..."
  let energyDiff = (energies.max - energies.min) / energies.len.float

  for R in 0 ..< nRadius:
    if R mod 10 == 0:
      echo &"Radius {R} of {nRadius}, {(R.float / nRadius.float) * 100:.2f} % done."
    n_eInt = n_es[R]
    radiiFlux = 0.0
    temperature = temperatures[R]
    let
      n_esR = n_es[R].float
      n_eFloat = n_esfloat[R] * 7.683e-24
      temp = temperatures[R].float
      radius = 0.0015 + R.float * 0.0005
      bfieldR = bfield(radius) #or something like that
      rho_keV = rho[R] * 7.683e-24 * 5.60958616722e29 #g/cm³ to keV⁴
    n_e_keV = pow(10.0, (n_esR * 0.25)) * 7.683e-24 # was 1/cm³ #correct conversion #now keV³
    n_e_keV = n_eFloat
    #echo n_e_keV, " ", n_eFloat
    let temp_keVTable = pow(10.0, (temp * 0.025)) * 8.617e-8 # was K # correct conversion
    let temp_keVFloat = tempFloat[R] * 8.617e-8
    temp_keV = temp_keVFloat
    var temp_K = pow(10.0, (temp * 0.025))
    let
      debye_scale_squared = (4.0 * PI * alpha / temp_keV) *
                              (n_e_keV + n_Z[R][1] * 7.645e-24 +
                               4.0 * n_Z[R][2] * 7.645e-24 )
      debye_scale = sqrt(debye_scale_squared)
      y = debye_scale / (sqrt( 2.0 * m_e_keV * temp_keV))
    for (iEindex, energy_keV) in pairs(energies):
      var
        sum = 0.0
        absCoef = 0.0
        table = 1.0
      let
        w = energy_keV / temp_keVTable #toFloat(dfMesh["u"][iE.int])

      #let n_bar_keV = n_Z[R][1] * 7.645e-24 + n_Z[R][2] * 7.645e-24 + alphaR[R] * metallicity(r)) * density(r)/((1.0E+9*eV2g)*atomic_mass_unit
      if w >= 20.0 or w <= 0.0732: #because the tables dont go beyond that, apparently because the axion production beyond that is irrelevant #except for He, maybe find a better solution
        for (Z_str, Z) in iterEnum(ElementKind):
          # TODO: avoid looping over unneeded elements here
          if Z in noElement:
            continue
          sum = sum + n_Z[R][Z] * 0.0
        absCoefs[R, iEindex] = sum * 1.97327e-8 * 0.528e-8 * 0.528e-8 *
                               (1.0 - exp(-energy_keV / temp_keV))
      else:
        table = spline.eval(w)
        let tempTab = opElements[temperature]
        for (Z_str, Z) in iterEnum(ElementKind):
          # TODO: avoid looping over unneeded elements here
          if Z in noElement:
            continue
          let opacity = tempTab[ElementKind(Z)].densityTab[n_eInt].interp.eval(table)
          # opacities in atomic unit for lenth squared: 0.528 x10-8cm * 0.528 x10-8cm = a0² # 1 m = 1/1.239841336215e-9 1/keV and a0 = 0.528 x10-10m
          if Z > 2:
            sum += n_Z[R][Z] * opacity

        absCoef = sum * 1.97327e-8 * 0.528e-8 * 0.528e-8 * (1.0 - exp(-energy_keV / temp_keV)) # is in keV

        absCoefs[R, iEindex] = absCoef
      ## Now it's left to calculate the emission rates
      ## making the same approximation as for n_e calculation

      ##ion density weighted by charge^2 from Raffelt
      let nZZ2_raffelt_keV = (rho[R] / amu) * 7.683e-24 #seems to be correct
      let ffterm = freefreeEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV, nZZ2_raffelt_keV, w, y)
      ## includes contribution from ff, fb and bb processes and a part of the Comption contribution ## keV³ / keV² = keV :

      let
        term1 = term1(g_ae, energy_keV, (absCoefs[R, iEindex]), e_charge, m_e_keV, temp_keV)
        term2 = term2(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV)# completes the Compton contribution #keV
        term3 = bremsEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV, w, y) # contribution from ee-bremsstahlung
        compton = comptonEmrate(alpha, g_ae, energy_keV, n_e_keV, m_e_keV, temp_keV)
      let
        primakoff = primakoff(temp_keV, energy_keV, gagamma, debye_scale_squared, alpha, n_e_keV, m_e_keV, n_Z[R][2], n_Z[R][1])
        longPlas = longPlasmon(energy_keV, n_e_keV, m_e_keV, alpha, bfieldR, temp_keV, (absCoefs[R, iEindex]), gagamma) #is correct with the absorbtion coefficient #https://arxiv.org/pdf/2006.10415.pdf make similar log pictures
        transPlas = transPlasmon(energy_keV, n_e_keV, m_e_keV, alpha, bfieldR, temp_keV, (absCoefs[R, iEindex]), gagamma) #is correct with the absorbtion coefficient
        iron57 = iron(ganuclei, temp_keV, energy_keV, rho_keV)
      term1s[R, iEindex] = term1
      comptons[R, iEindex] = compton
      term3s[R, iEindex] = term3
      ffterms[R, iEindex] = ffterm
      primakoffs[R, iEindex] = primakoff
      longPlasmons[R, iEindex] = longPlas
      transPlasmons[R, iEindex] = transPlas
      iron57s[R, iEindex] = iron57
      abc[R, iEindex] = compton +  term1 + term3 + ffterm
      let total_emrate = compton +  term1 + term3 + ffterm  + transPlas + primakoff  + longPlas + iron57
      let total_emrate_s = total_emrate / (6.58e-19) # in 1/sec
      emratesS[R, iEindex] = total_emrate
      emratesInS[R, iEindex] = total_emrate_s

      radiiFlux += totalEmRates * energyDiff #iron57 * 0.001
      #if w <= 0.0732 :
        #echo transPlas
      # if want to have absorbtion coefficient of a radius and energy: R = (r (in % of sunR) - 0.0015) / 0.0005


    totalRadiiFlux[R] = radiiFlux
    radii[R] = radius
    #echo R, " ", omegaPlasmonSq(alpha, n_e_keV, m_e_keV)

  #echo longPlasmons[0]
  echo "creating all plots..."

  #[var dfNZ = seqsToDf({ "Radius": rs,
                        "nZ": nZs,
                        "Z": zs})
  echo dfNZ
  ggplot(dfNZ, aes("Radius", "nZ", color = "Z")) +
    geom_point() +
    ggtitle("Radius versus atomic density for different Z") +
    ggsave("out/radius_nZ_Z.pdf")]#

  #var dfOp = seqsToDf({ "Radius": rs2,
  #                      "opacity": ops,
  #                      "energies": engs})
  #echo dfOp
  #ggplot(dfOp, aes("Radius", "opacity", color = "energies")) +
  #  geom_point() +
  #  ggtitle("Radius versus opacity for different energies") +
  #  ggsave("out/radius_op_energies.pdf")
  let dfRadii = seqsToDf({"Radius": radii,
                        "Flux": totalRadiiFlux})
  echo dfRadii

  ggplot(dfRadii, aes("Radius", "Flux")) +
    geom_line() + #size = some(0.5)
    xlab("Solar radius") +
    ylab("Flux") +
    ggtitle(&"Differential solar axion flux for g_ae = {g_ae}, g_aγ = {g_agamma} GeV⁻¹, g_aN = {ganuclei}") +
    ggsave("out/radFlux.pdf", width = 800, height = 480)

  var diffFluxDf = newDataFrame()
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, emratesS, "Total flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, term1s, "FB BB Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, comptons, "Compton Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, term3s, "EE Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, ffterms, "FF Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, primakoffs, "Primakoff Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, longPlasmons, "LP Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, transPlasmons, "TP Flux")
  diffFluxDf.add getFluxFractionR(energies, df, n_es, temperatures, iron57s, "57Fe Flux")
  echo diffFluxDf
  let
    ironflux = getFluxFractionR(energies, df, n_es, temperatures, iron57s, "57Fe Flux")
  echo ironflux
  let
    ir = ironflux["diffFlux"].toTensor(float).torawseq
  echo "57Fe Flux ", ir.foldl(a + b) * 0.001 / ganuclei / ganuclei * 3.171e-12 , " g_aN² cm⁻2 s⁻1"
  #let diffFluxDf = seqsToDf({ "Energy / eV" : energies,
  #                            "Flux / keV⁻¹ m⁻² yr⁻¹" : diff_fluxs })
  #diffFluxDf.write_csv(&"axion_diff_flux_gae_{g_ae}_gagamma_{g_agamma}.csv")
  #let totalFlux = simpson(diff_fluxs, energies.mapIt(it * 1e-3))
  #echo "The total axion Flux in 1/(y m^2):", totalFlux

  #[let dfEmrate = seqsToDf({ "energy": energies,
                            "emrate": emratesS[4, _].squeeze.clone })
    .mutate(f{"flux" ~ `emrate` * `energy` * `energy` * 0.5 / Pi / Pi})
  ggplot(dfEmrate, aes("energy", "flux")) +
    geom_line() +
    ggsave("out/emrate_R10.pdf")

  let dfAbscoef = seqsToDf({ "energy": energies,
                            "absCoefs": absCoefs[10, _].squeeze.clone })
  ggplot(dfAbscoef, aes("energy", "absCoefs")) +
    geom_line() +
    ggsave("out/abscoefs_R10.pdf")


  ggplot(diffFluxDf, aes("Energy", "Flux", color = "Radius")) +
    geom_line() +
    xlab("Axion energy [keV]") +
    ylab("Flux [keV⁻¹ y⁻¹ m⁻²]") +
    ggtitle(&"Differential solar axion flux for g_ae = {g_ae}, g_aγ = {g_agamma} GeV⁻¹") +
    margin(right = 6.5) +
    ggsave("out/diffFlux_radii.pdf", width = 800, height = 480)

  when false:
    let dfDiffflux = seqsToDf({ "Axion energy [eV]": energieslong,
                                "Fluxfraction [keV⁻¹y⁻¹m⁻²]": fluxes,
                                "type": kinds })]#

  ggplot(difffluxDf, aes("Energy", "diffFlux", color = "type")) +
    geom_line() + #size = some(0.5)
    xlab("Axion energy [eV]") +
    ylab("Flux [keV⁻¹ y⁻¹ m⁻²]") +
    #ylim(0, 2.5e24) +
    #xlim(0.0, 1000.0) +
    xlim(0.0, 15.0) +
    #scale_y_log10() + #
    scale_y_continuous() +
    #scale_x_log10() +
    ggtitle(&"Differential solar axion flux for g_ae = {g_ae}, g_aγ = {g_agamma} GeV⁻¹, g_aN = {ganuclei}") +
    margin(right = 6.5) +
    ggsave("out/diffFlux.pdf", width = 800, height = 480)

  result = emratesS

when isMainModule:

  let solarModel = main()
  echo "writing to csv"
  solarModel.to_csv("solar_model_tensor15.csv")

  when false:
    ## TODO: better approach to store the full "solar model" as a CSV from a DF
    var dfSolarModel = newDataFrame()
    var radii = newSeq[string]()
    for R in 0 ..< solarModel.len:
      let r = $R
      dfSolarModel[r] = solarModel[R]
      radii.add r
    dfSolarModel["Energy"] = energies

    dfSolarModel = dfSolarModel.gather(radii, key = "Radii", value = "EmissionRates")
    # energy / float, Radius / string!, emission rate / float
    dfSolarModel = dfSolarModel.mutate(f{string -> int: "Radii" ~ parseInt(df["Radii"][idx])})
      .mutate(f{int -> float: "Radii" ~ `Radii`.float * 0.0005 + 0.0015})
    # energy / float, Radius / float, emission rate / float
    echo dfSolarModel
    dfSolarModel.write_csv("created_solar_model_em.csv")
