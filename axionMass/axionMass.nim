import sequtils, seqmath, ggplotnim, strformat, algorithm, nlopt, options, strutils

proc effPhotonMass*(ne: float): float =
  ## returns the effective photon mass for a given electron number density
  const alpha = 1.0 / 137.0
  const me = 511e3 # 511 keV
  # note the 1.97e-7 cubed to account for the length scale in `ne`
  result = sqrt( pow(1.97e-7, 3) * 4 * PI * alpha * ne / me )

proc numDensity(c: float): float =
  ## converts a molar concentration in mol / m^3 to a number density
  const Z = 2 # number of electron in helium atom
  result = Z * 6.022e23 * c

proc molarAmount(p, vol, temp: float): float =
  ## calculates the molar amount of gas given a certain pressure,
  ## volume and temperature
  ## the pressure is assumed in mbar
  const gasConstant = 8.314 # joule K^-1 mol^-1
  let pressure = p * 1e2 # pressure in Pa
  result = pressure * vol / (gasConstant * temp)
  #echo "Molar amount for P = ", pressure, " Pa is ", result

proc babyIaxoEffMass(p: float): float =
  ## calculates the effective photon (and thus axion mass) for BabyIAXO given
  ## a certain helium pressure in the BabyIAXO magnet
  const vol = 10.0 * (PI * pow(0.3, 2)) # 10m length, bore radius 30 cm
  # UPDATE: IAXO will be run at 4.2 K instead of 1.7 K
  # const temp = 1.7 # assume 1.7 K same as CAST
  const temp = 4.2
  once:
    echo "BabyIAXO magnet volume is ", vol, " m^3"
    echo "BabyIAXO magnet temperature is ", temp, " K"
  let amountMol = molarAmount(p, vol, temp) # amount of gas in mol
  let numPerMol = numDensity(amountMol / vol) # number of electrons per m^3
  #echo "Num electrons per m^3 ", numPerMol
  result = effPhotonMass(numPerMol)

proc effPhotonMass2(p : float, length : float, radBore : float, temp : float): float =
  ## calculates the effective photon (and thus axion mass) given
  ## a certain helium pressure in the magnet
  ## length and radBore in m
  let vol = length * (PI * pow(radBore, 2)) # 10m length, bore radius 30 cm
  # UPDATE: IAXO will be run at 4.2 K instead of 1.7 K
  # const temp = 1.7 # assume 1.7 K same as CAST
  let amountMol = molarAmount(p, vol, temp) # amount of gas in mol
  let numPerMol = numDensity(amountMol / vol) # number of electrons per m^3
  #echo "Num electrons per m^3 ", numPerMol
  result = effPhotonMass(numPerMol)

proc logspace(start, stop: float, num: int, base = 10.0): seq[float] =
  ## generates evenly spaced points between start and stop in log space
  let linear = linspace(start, stop, num)
  result = linear.mapIt(pow(base, it))

proc makePlot(pstart, pstop: float, fname: string, log = false) =
  let pressures = logspace(pstart, pstop, 1000) # 1000 values from 1e-5 mbar to 500 mbar
  let masses = pressures.mapIt(babyIaxoEffMass(it)) # corresponding masses
  # convert both seqs to a dataframe
  let df = seqsToDf({"P / mbar" : pressures, "m_a / eV" : masses})
  let plt = ggplot(df, aes("P / mbar", "m_a / eV")) +
    geom_line() +
    ggtitle("Sensitive axion mass in eV depending on helium pressure in mbar")
  if not log:
    plt + ggsave(fname)
  else:
    plt + scale_x_log10() + scale_y_log10() + ggsave(fname)

#makePlot(-6.0, 2.0, "axionMassesFullRange.pdf")
#
#makePlot(-6.0, -2.0, "axionMassesZoomed.pdf")
#
#makePlot(-6.0, 2.0, "axionMassesLogLog.pdf", log = true)

proc analyticalCalc(p: float): float =
  ## calculates the effective mass for a given pressure from a single equation
  ## derived from inserting all above into one
  result = 1.94081e-2 * sqrt(4 * PI * p )

#block:
#  let pressures = logspace(-6.0, 2.0, 1000)
#  let massesNum = pressures.mapIt(babyIaxoEffMass(it))
#  let massesAna = pressures.mapIt(analyticalCalc(it))
#  #echo massesNum
#  #echo massesAna
#  func almostEqual(a, b: seq[float]): bool = zip(a, b).mapIt(abs(it[0] - it[1]) < 1e-5).allIt(it)
#  doAssert massesNum.almostEqual(massesAna)

proc vacuumMassLimit(energy: float, magnetLength = 10.0): float =
  ## given an axion energy in keV, calculate the limit of coherence
  ## for the vacuum case in BabyIAXO
  # note the length scale 1.97e-7 to take into account the meter scale in
  # babyIaxoLen
  let babyIaxoLen = magnetLength / 1.97e-7 # length in "eV"
  result = sqrt(PI * 2 * energy * 1e3 / babyIaxoLen) # convert energy to eV

#let energies = linspace(0.0, 10.0 , 1000) # from 0 to 10 keV
#let massLimits = energies.mapIt(vacuumMassLimit(it))
#let df = seqsToDf({"E / keV" : energies, "m_a / eV" : massLimits})
#ggplot(df, aes("E / keV", "m_a / eV")) +
#  geom_line() +
#  ggtitle("Sensitive axion mass limit in eV for BabyIAXO in vacuum run") +
#  ggsave("vacuumMassLimit.pdf")

#let iaxoLimit = vacuumMassLimit(4.2, 20.0)
#func almostEqual(a, b: float, eps = 1e-3): bool = abs(a - b) < eps
#doAssert iaxoLimit.almostEqual(1.61e-2) # value from IAXO gas phase study
#echo "Full IAXO mass limit @ 4.2 keV = ", iaxoLimit, " eV"
#
#const babyIaxoVacuumMassLimit = vacuumMassLimit(4.2)
#echo "BabyIAXO mass limit @ 4.2 keV = ", babyIaxoVacuumMassLimit, " eV"

proc logMassAttenuation(e: float): float =
  ## calculates the logarithm of the mass attenuation coefficient for a given
  ## energy `e` in `keV` and the result in `cm^2/g`
  result = -1.5832 + 5.9195 * exp(-0.353808 * e) + 4.03598 * exp(-0.970557 * e)

#let logMuOverRho = energies.mapIt(logMassAttenuation(it))
## now the non-log values
#let muOverRho = logMuOverRho.mapIt(exp(it))
#[
const massAttenuationFile = "mass_attenuation_nist_data.txt"
# skip one line after header, second header line
var dfMuRhoTab = toDf(readCsv(massAttenuationFile, skipLines = 1, sep = ' ', header = "#"))#
  # convert MeV energy to keV
  .mutate(f{"Energy" ~ "Energy" * 1000.0})
  .filter(f{"Energy" >= energies.min and "Energy" <= energies.max})

# create df of interpolated values
let dfMuRhoInterp = seqsToDf({ "E / keV" : energies,
                               "mu/rho" : muOverRho,
                               "log(mu/rho)" : logMuOverRho})
# rename the columns of the tabulated values df and create a log column
dfMuRhoTab = dfMuRhoTab.rename(f{"E / keV" ~ "Energy"})
  .mutate(f{"log(mu/rho)" ~ ln("mu/rho")})
# build combined DF
let dfMuRho = bind_rows([("Interpolation", dfMuRhoInterp),
                         ("NIST", dfMuRhoTab)],
                        id = "type")
echo dfMuRho
# create plot of log values
ggplot(dfMuRho, aes("E / keV", "log(mu/rho)", color = "type")) +
  geom_line(data = dfMuRho.filter(f{"type" == "Interpolation"})) +
  geom_point(data = dfMuRho.filter(f{"type" == "NIST"})) +
  ggtitle("Mass attenuation coefficient interpolation and data") +
  ggsave("log_mass_attenuation_function.pdf")

# and the plot of the raw mu/rho values
ggplot(dfMuRho, aes("E / keV", "mu/rho", color = "type")) +
  geom_line(data = dfMuRho.filter(f{"type" == "Interpolation"})) +
  geom_point(data = dfMuRho.filter(f{"type" == "NIST"})) +
  scale_y_log10() +
  ggtitle("Mass attenuation coefficient interpolation and data") +
  ggsave("mass_attenuation_function.pdf")
]#
proc momentumTransfer(m_gamma, m_a: float, E_a = 4.2): float =
  ## calculates the momentum transfer for a given effective photon
  ## mass `m_gamma` and axion mass `m_a` at an axion energy of
  ## 4.2 keV `E_a` (by default).
  #const c = 299792458
  result = abs((m_gamma * m_gamma - m_a * m_a) / (2 * E_a * 1000.0))

proc densityTemp(p: float, temp : float): float =
  ## returns the density of the gas for the given pressure.
  ## The pressure is assumed in `mbar` and the temperature (in `K`).
  ## The default temperature corresponds to BabyIAXO aim.
  ## Returns the density in `g / cm^3`
  const gasConstant = 8.314 # joule K^-1 mol^-1
  const M = 4.002602 # g / mol
  let pressure = p * 1e2 # pressure in Pa
  # factor 1000 for conversion of M in g / mol to kg / mol
  result = pressure * M / (gasConstant * temp * 1000.0)
  # convert the result to g / cm^3 for use with mass attenuations
  result = result / 1000.0

proc axionConversionProb(m_a, m_gamma, gamma: float, length = 10.0): float =
  ## given an axion mass and an inverse absorption length
  ## calculates the conversion probabilty with BabyIAXOs magnet
  ## properties. Axion coupling constant taken to be `1` +1e-11+.
  ## By default uses BabyIAXO length of `10m`
  # both `g_agamma` and `B` only scale the absolute value `P`, does not matter
  const g_agamma = 1.0 #1e-12
  const B = 4.5 # T, actually don't know the real number right now #it's 5.4 T, I think
  # convert length in `m` to `eV`
  let L = length / 1.97e-7 # m
  let q = momentumTransfer(m_gamma, m_a)
  let term1 = pow(g_agamma * B / 2.0, 2)
  let term2 = 1.0 / (q * q + gamma * gamma / 4)
  let term3 = 1.0 + exp(-gamma * L) - 2 * exp(-gamma * L / 2) * cos(q * L)
  result = term1 * term2 * term3

proc axionConversionProb2*(m_a : float, energyAx : float, pressure : float, temp : float, length : float, radBore : float, g_agamma : float, B : float): float =
  ## given an axion mass and an inverse absorption length
  ## calculates the conversion probabilty
  # both `g_agamma` and `B` only scale the absolute value `P`, does not matter
  let energyGamma = energyAx #in keV
  ## for a given pressure in `mbar` returns the attenuation length
  ## for a given pressure in `mbar` returns the attenuation length
  ## `Γ` in units of `eV`.
  # multiply by `100` to convert `Γ` from `1 / cm` to `1 / m`
  # multiply by `1.97e-7` to convert `1 / m` to `1 / eV`
  # Temperature in K
  let gamma = 1.97e-7 * 100.0 * densityTemp(pressure , temp) * exp(logMassAttenuation(energyGamma))
  let m_gamma = effPhotonMass2(pressure, length, radBore, temp)
  #echo m_gamma
  #echo m_gamma - m_a
  #echo "m_a ", m_a
  #echo "gamma ", gamma
  # convert length in `m` to `1/eV`
  let L = length / 1.97e-7
  let g_agammaEV = g_agamma * 1e-9 #usually in units of 1/GeV but need units of 1/eV
  # convert B from T to eV^2
  let beV = B * 1e3 / 1.444
  let q = momentumTransfer(m_gamma, m_a, energyAx)
  let term1 = pow(g_agammaEV * beV / 2.0, 2)
  let term2 = 1.0 / (q * q + gamma * gamma / 4)
  let term3 = 1.0 + exp(-gamma * L) - 2 * exp(-gamma * L / 2) * cos(q * L)
  result = term1 * term2 * term3

#let pressures = logspace(-6.0, 2.0, 4) # take only few values for a start
#let m_gammas = pressures.mapIt(babyIaxoEffMass(it))

proc density(p: float, temp = 4.2): float =
  ## returns the density of the gas for the given pressure.
  ## The pressure is assumed in `mbar` and the temperature (in `K`).
  ## The default temperature corresponds to BabyIAXO aim.
  ## Returns the density in `g / cm^3`
  const gasConstant = 8.314 # joule K^-1 mol^-1
  const M = 4.002602 # g / mol
  let pressure = p * 1e2 # pressure in Pa
  # factor 1000 for conversion of M in g / mol to kg / mol
  result = pressure * M / (gasConstant * temp * 1000.0)
  # convert the result to g / cm^3 for use with mass attenuations
  result = result / 1000.0


proc attenuationLength(p: float): float =
  ## for a given pressure in `mbar` returns the attenuation length
  ## `Γ` in units of `eV`.
  # multiply by `100` to convert `Γ` from `1 / cm` to `1 / m`
  # multiply by `1.97e-7` to convert `1 / m` to `1 / eV`
  result = 1.97e-7 * 100.0 * density(p) * exp(logMassAttenuation(4.2))

#let densities = pressures.mapIt(density(it))
#
#let gammas = densities.mapIt(it * exp(logMassAttenuation(4.2)))

proc genAxionConvPlot(gamma, m_gamma: float, nameSuffix = "",
                      start = 0.0, stop = 0.0, length = 10.0) =
  ## generates a single axion conversion probabilit plot
  # calculate reasonable `m_a` values
  echo "Gamma: ", gamma
  echo "m_gamma: ", m_gamma
  var m_a: seq[float]
  if start != stop:
    m_a = linspace(start, stop, 1000)
  else:
    m_a = linspace(m_gamma * 0.99, m_gamma * 1.01, 1000)

  let qs = m_a.mapIt(momentumTransfer(m_gamma, it))
  # plot momentum transfers by themselves
  let dfMom = seqsToDf(m_a, qs)
  ggplot(dfMom, aes("m_a", "qs")) +
    geom_line() +
    ggsave("momentumTransfers.pdf")

  let prob = m_a.mapIt(axionConversionProb(it, m_gamma, gamma, length = length))
  let df = seqsToDf({ "m_a / eV" : m_a,
                      "P_a->gamma" : prob })
  echo df
  #echo df.summarize(f{"P_a->gamma" ~ min("P_a->gamma")})
  #echo df.pretty(-1)
  ggplot(df, aes("m_a / eV", "P_a->gamma")) +
    geom_line() +
    ggtitle(&"Axion conversion probability for Γ = {gamma:.2e}, m_γ = {m_gamma:.2e}") +
    ggsave(&"axion_conversion_prob_{nameSuffix}.pdf")
  if nameSuffix == "1":
    writeFile("dfdata.txt", df.pretty(-1))

#doAssert gammas.len == m_gammas.len
#for i in 0 ..< gammas.len:
#  genAxionConvPlot(gammas[i], m_gammas[i], $i)
#let dftest = seqsToDf({ "gammas" : gammas,
#                        "m_gammas" : m_gammas })
#ggplot(dftest, aes("gammas", "m_gammas")) +
#  geom_line() +
#  ggsave("test_gamma.pdf")

proc pressureAtTemp(p: float, T = 4.2, T2 = 293.0): float =
  ## converts a given pressure `p` at `T` to an equivalent pressure at `T2`
  result = p * T2 / T

#let p1bar = pressureAtTemp(1000.0, 293, 4.2)
#let p3bar = pressureAtTemp(3000.0, 293, 4.2)
#echo "Pressure for equivalent of 1 bar @ 293 K ", p1bar
#echo "Pressure for equivalent of 3 bar @ 293 K ", p3bar
#let m_gamma_1bar = babyIaxoEffMass(p1bar)
#let m_gamma_3bar = babyIaxoEffMass(p3bar)
#echo "m_gamma for pressure equivalent of 1 bar @ 293 K ", m_gamma_1_bar
#echo "m_gamma for pressure equivalent of 3 bar @ 293 K ", m_gamma_3_bar
#doAssert almostEqual(m_gamma_1_bar, 0.26)
#doAssert almostEqual(m_gamma_3_bar, 0.4483, eps = 1e-2)
#
#let gamma_1bar = attenuationLength(p1bar)
#let gamma_3bar = attenuationLength(p3bar)
#genAxionConvPlot(gamma_1bar, m_gamma_1bar, "1bar_equiv")#, 0.2593, 0.2616
#genAxionConvPlot(gamma_3bar, m_gamma_3bar, "3bar_equiv")#, 0.4485, 0.4535
## and now for reference the IAXO (`20m`) plots
#genAxionConvPlot(gamma_1bar, m_gamma_1bar, "1bar_equiv_20m", length = 20.0)
#genAxionConvPlot(gamma_3bar, m_gamma_3bar, "3bar_equiv_20m", length = 20.0)
#
import mpfit
proc gaussFit(p_ar: seq[float], x: float): float =
  result = p_ar[0] * gauss(x = x, mean = p_ar[1], sigma = p_ar[2])

proc calcConvProbCurve(gamma, m_gamma, pressure: float, length = 10.0,
                       massRange = none[tuple[low, high: float]]()):
    tuple[m_a, prob: seq[float]] =
  var m_a: seq[float]
  if massRange.isNone:
    if pressure < 0.01:
      echo "PRressure ", pressure
      m_a = linspace(m_gamma * 0.2, m_gamma * 1.8, 1000)
    else:
      m_a = linspace(m_gamma * 0.5, m_gamma * 1.5, 1000)
  else:
    m_a = linspace(massRange.get.low, massRange.get.high, 1000)
  let qs = m_a.mapIt(momentumTransfer(m_gamma, it))
  let prob = m_a.mapIt(axionConversionProb(it, m_gamma, gamma, length = length))
  result = (m_a: m_a, prob: prob)

proc fitToConvProb(gamma, m_gamma, pressure: float, nameSuffix = "",
                   length = 10.0, createPlot = true): (float, seq[float]) =
  ## generates a single axion conversion probabilit plot
  let (m_a, prob) = calcConvProbCurve(gamma, m_gamma, pressure, length)
  let p0 = @[prob.max, m_a[prob.argmax], 0.02]
  let (pRes, res) = fit(gaussFit, p0, m_a, prob, prob.mapIt(1.0))
  #echoResult(pRes, res = res)
  let probFit = m_a.mapIt(gaussFit(pRes, it))
  if createPlot:
    let df = seqsToDf({ "m_a / eV" : m_a,
                        "P_a->gamma" : prob,
                        "P_fit" : probFit })
    #echo df
    ggplot(df, aes("m_a / eV", "P_a->gamma")) +
      geom_line() +
      geom_line(aes(y = "P_fit", color = "Gaussian fit")) +
      ggtitle(&"Gaussian fit to P_a->gamma at p = {pressure:.4f} mbar") +
      ggsave(&"conv_prob_fit_{nameSuffix}.pdf")
  result = (pressure, pRes)

#let (drop, pRes1bar) = fitToConvProb(gamma_1bar, m_gamma_1bar, p1bar, "1bar")

proc fwhm(sigma: float): float =
  ## returns the FWHM of a gaussian for a `sigma`
  result = 2 * sqrt(2 * ln(2.0)) * abs(sigma)
#echo "FWHM @ 1bar room temperature = ", fwhm(pRes1bar[2])

#let pressuresFine = logspace(-6.0, 2.0, 1000)
#let gammasFine = pressuresFine.mapIt(attenuationLength(it))
#let mgammasFine = pressuresFine.mapIt(babyIaxoEffMass(it))
#var fwhmFine: seq[float]
when false:
  # NOTE: set to `true` if you want to recreate the fits and plots!
  for i in 0 ..< pressuresFine.len:
    let (p, pResI) = fitToConvProb(gammasFine[i], mgammasFine[i],
                                   pressuresFine[i], &"{pressuresFine[i]:.6f}",
                                   createPlot = true)
    fwhmFine.add fwhm(pResI[2])

when false:
  # NOTE: set to `true` if you want to redo these calculations
  var dfFwhm = seqsToDf({ "Pressure / mbar" : pressuresFine,
                          "FWHM / eV" : fwhmFine })
  dfFwhm = dfFwhm.mutate(f{"FWHM / eV" ~ abs("FWHM / eV")})
  ggplot(dfFwhm, aes("Pressure / mbar", "FWHM / eV")) +
    geom_point() +
    ylab(margin = 1.5) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle("FWHM / eV of a->gamma probability depending on pressure / mbar") +
    ggsave("fwhm_vs_pressure_full_loglog.pdf")
  # and a subset filtered to above 0.01 mbar in linear
  dfFwhm = dfFwhm.filter(f{"Pressure / mbar" > 0.01})
  ggplot(dfFwhm, aes("Pressure / mbar", "FWHM / eV")) +
    geom_point() +
    scale_y_log10() +
    ggtitle("FWHM / eV of a->gamma probability depending on pressure / mbar") +
    ggsave("fwhm_vs_pressure_ylog.pdf")

  let elements = pressuresFine.len
  doAssert pressuresFine.len == elements
  doAssert gammasFine.len == elements, " was " & $gammas.len
  doAssert mgammasFine.len == elements, " was " & $mgammas.len
  doAssert fwhmFine.len == elements, " was " & $fwhmFine.len

# finally write the data to a file
proc writeCsv(pressures, gammas, mgammas, fwhm: seq[float])
when false:
  writeCsv(pressuresFine, gammasFine, mgammasFine, fwhmFine)

proc writeCsv(pressures, gammas, mgammas, fwhm: seq[float]) =
  var f = open("fwhm_results.csv", fmWrite)
  f.write("# Pressures / mbar\t Γ / eV\t m_γ / eV & fwhm / eV\n")
  let elements = pressures.len
  var line = ""
  for i in 0 ..< elements:
    line = &"{pressures[i]},{gammas[i]},{mgammas[i]},{abs(fwhm[i])}\n"
    f.write(line)
  f.close()

#echo "Room temp: ", pressures.mapIt(pressureAtTemp(it) / 1000.0)

proc m_gamma_alternative(p: float): float =
  # gallard alternative m_gamma. Says 0.02, more exact is:
  result = sqrt(0.01988 * p / 4.2)

#let compPress = pressureAtTemp(1000, 293.0, 4.2)
#echo "cryo pressure for 1 bar at room: ", compPress
## get m_a at compPress
#echo "m_a @ 1 bar @ room temp: ", babyIaxoEffMass(compPress)
#echo "m_a alt. @ 1 bar room temp ", m_gamma_alternative(compPress)
#doAssert babyIaxoEffMass(compPress).almostEqual(m_gamma_alternative(compPress))


proc intensitySuppression*(energy, distance, pressure: float, temp = 293.15): float =
  ## calculates the suppression factor of the intensity of X-rays
  ## with `energy` (in keV) after `distance` (in `m`) at a helium `pressure`
  ## (in `mbar`) and temperature of the gas of `temp` (in `K`).
  ## By default assume room temperature for calc of beam line filled with
  ## gas.
  let massAtt = exp(logMassAttenuation(energy))
  # calculate density from pressure
  let rho = density(pressure, temp)
  # factor 100 is to convert `distance` to `cm` from `m`
  result = exp(-massAtt * rho * distance * 100)

proc intensitySuppression2*(energy : float, distanceMagnet : float, distancePipe : float, pressure: float, tempMagnet : float, tempPipe : float64): float =
  ## calculates the suppression factor of the intensity of X-rays
  ## with `energy` (in keV) after `distance` (in `m`) at a helium `pressure`
  ## (in `mbar`) and temperature of the gas of `temp` (in `K`).
  ## By default assume room temperature for calc of beam line filled with
  ## gas. temp = 293.15
  let massAtt = exp(logMassAttenuation(energy))
  # calculate density from pressure
  let rhoMagnet = density(pressure, tempMagnet)
  let rhoPipe = density(pressure, tempPipe)
  # factor 100 is to convert `distance` to `cm` from `m`
  result = exp(-massAtt * rhoPipe * distancePipe * 100) * exp(-massAtt * rhoMagnet * distanceMagnet * 100)

proc calcFilledBeamline(length: float, energies: seq[float],
                        toPlot = true): DataFrame =
  ## calculates the X-ray suppression of a beamline filled with gas at
  ## `1 bar`, `293.15 K` of `5 m` (default) length
  # take intensity suppression as is, assumes incoming intensity is `1`
  let p = 1000.0 # mbar
  let intensities = energies.mapIt(intensitySuppression(it, length, p))
  result = seqsToDf({ "E / keV" : energies,
                      "Transmission" : intensities })
  if toPlot:
    ggplot(result, aes("E / keV", "Transmission")) +
      geom_line() +
      ggtitle(&"Transmission after {length} m of He at 1 bar, 293.15 K") +
      ggsave(&"transmission_beamline_he_{length}_m.pdf")

#let df5m = calcFilledBeamline(5.0, energies)
#let df7_5m = calcFilledBeamline(7.5, energies)
#echo "\n\n\n"
#echo df5m.tail(10)
#echo df7_5m.tail(10)
#echo "done\n\n\n"
#[
proc transmissionsPlusSpectrum =
  ## creates a plot of the expected transmissions for filled beamlines,
  ## the polypropylene window and the axion electron spectrum.
  let polypropDf = toDf(readCsv("polypropylene_window_10micron.txt", sep = ' ', skipLines = 1,
                                header = "#"))
    .mutate(f{"E / keV" ~ "Energy/eV" / 1000.0})
    .rename(f{"10 µm PP" ~ "Transmission"})
  # now using the `polyPropDf` energies, calculate the data frames for the He losses, so that
  # we have the exact same energies and we can easily multiply them afterwards
  let energiesPP = polyPropDf["E / keV"].vToSeq.mapIt(it.toFloat)
  let df5m = calcFilledBeamline(5.0, energiesPP).rename(f{"5 m He" ~ "Transmission"})
  let df7_5m = calcFilledBeamline(7.5, energiesPP).rename(f{"7.5 m He" ~ "Transmission"})
  # do an ugly multi join
  var dfPPHe = inner_join(polyPropDf, df5m, by = "E / keV").innerJoin(df7_5m, by = "E / keV")
  # finally add polyProp * He
  dfPPHe = dfPPHe.mutate(f{"5mHe+PP" ~ "10 µm PP" * "5 m He"},
                         f{"7.5mHe+PP" ~ "10 µm PP" * "7.5 m He"})
    .gather(["10 µm PP", "5 m He", "5mHe+PP", "7.5 m He", "7.5mHe+PP"],
            key = "Setup", value = "Transmission")
  let gaeDf = toDf(readCsv("axion_gae_flux.dat", sep = ' ', skipLines = 9, header = "#"))
    .rename(f{"E / keV" ~ "Energy/keV"})
    .mutate(f{"Phi" ~ "Phi" / max("Phi") * 100.0})
    .filter(f{"E / keV" <= 10.0})
  var fullDf = bind_rows([("Axion-Electron flux", gaeDf),
                          ("Setups", dfPPHe)],
                         id = "Type")
    .select("Type", "E / keV", "Phi", "Transmission", "Setup")
    .mutate(f{"Transmission" ~ "Transmission" * 100.0})
  echo fullDf.pretty(-1)
  echo fullDf
  ggplot(fullDf.filter(f{"Type" != "Axion-Electron flux"}), aes(x = "E / keV", y = "Transmission")) +
    geom_line(aes(color = "Setup")) +
    geom_line(aes(x = "E / keV", y = "Phi"),
              data = fullDf.filter(f{"Type" == "Axion-Electron flux"})) +
    scale_y_continuous("Transmission / %", secAxis = secAxis(name = "Axion flux / a.u.")) +
    legend_position(x = 0.88, y = 0.1) +
    #xlim(0.0, 3.0, outsideRange = "drop") +
    ggtitle("Comparison of He filled beamline (1 bar, 293 K) and 10 µm window") +
    ggsave("window_he_transmissions_axion_flux.pdf", width = 853, height = 480)

transmissionsPlusSpectrum()
]#
proc vacuumConversionProb(m_a: float, length = 10.0): float =
  ## calculates the conversion probability in BabyIAXO for the given axion
  ## mass `m_a`
  # both `g_agamma` and `B` only scale the absolute value `P`, does not matter
  const g_agamma = 1.0 #1e-11
  const B = 4.5 # T, actually don't know the real number right now
  # convert length in `m` to `eV`
  let L = length / 1.97e-7 # m
  let E_a = 4.2 # keV mean value for Primakoff spectrum
  let q = momentumTransfer(m_gamma = 0.0, m_a = m_a)
  let term1 = pow(g_agamma * B * L / 2.0, 2)
  let term2 = pow(sin(q * L / 2.0) / (q * L / 2.0), 2.0)
  result = term1 * term2

proc vacuumConversionProbApprox(m_a: float, length = 10.0): float =
  ## approximates the vacuum conversion probability to second order of taylor
  ## expansion (4th power of `sin^2` argument)
  # both `g_agamma` and `B` only scale the absolute value `P`, does not matter
  const g_agamma = 1.0 #1e-11
  const B = 4.5 # T, actually don't know the real number right now
  # convert length in `m` to `eV`
  let L = length / 1.97e-7 # m
  let E_a = 4.2 # keV mean value for Primakoff spectrum
  let q = momentumTransfer(m_gamma = 0.0, m_a = m_a)
  let term1 = pow(g_agamma * B * L / 2.0, 2)
  let term2 = 1.0 - pow(q * L / 2.0, 2.0) / 3.0
  result = term1 * term2

proc vacuumConvProbPlot(useApprox = false) =
  ## generates the typical plot of the conversion probability in the range
  ## 1 µeV to 1 eV for g_a,gamma = 1.0 (for simplicity)
  # we use linspace so we have more fidelity where the curve is more interesting
  let vacMasses = linspace(1e-6, 1.0, 1000)
  var probs: seq[float]
  var outname = ""
  if not useApprox:
    probs = vacMasses.mapIt(vacuumConversionProb(it))
    outname = "vacuum_conversion_prob.pdf"
  else:
    probs = vacMasses.mapIt(vacuumConversionProbApprox(it))
    outname = "vacuum_conversion_prob_taylor.pdf"
  let dfVac = seqsToDf({ "m_a / eV" : vacMasses,
                         "P_a,gamma" : probs})

  ggplot(dfVac, aes("m_a / eV", "P_a,gamma")) +
    geom_line() +
    #geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    ggsave(outname)
#vacuumConvProbPlot()
# The below results in negative values for the probability, which breaks the
# log log plot. Why negative from approximation?
#vacuumConvProbPlot(useApprox = true)

proc pressureToGammaMGamma(pressure: float): tuple[gamma, m_gamma: float] =
  ## returns a tuple of the corresponding `gamma` and
  ## `m_gamma` values associated with the pressure for BabyIAXO filled with
  ## He.
  result = (gamma: attenuationLength(pressure),
            m_gamma: babyIaxoEffMass(pressure))

proc getDfVac(massRange: Option[tuple[low, high: float]]): DataFrame =
  let r = if massRange.isNone: (low: 1e-6, high: 1.0) else: massRange.get
  let vacMasses = linspace(r.low, r.high, 1000)
  let vacProbs = vacMasses.mapIt(vacuumConversionProb(it))
  result = seqsToDf({ "m_a / eV" : vacMasses,
                      "P_a->gamma" : vacProbs})
  # echo result.pretty(50)

type
  MassFit = object
    discard
proc vacuumMassGivenProb(prob: float): float =
  ## determines the mass corresponding to the given probability.
  ## Since `P` even for the vacuum case is not really invertible,
  ## we'll do it numerically.
  proc probDiff(p0: seq[float], dummy: MassFit): float =
    result = abs(vacuumConversionProb(p0[0]) - prob)
  var opt = newNloptOpt[MassFit](LN_COBYLA, 1)
  # hand the function to fit as well as the data object we need in it
  var varStr = newVarStruct(probDiff, MassFit())
  opt.setFunction(varStr)
  let (params, minVal) = opt.optimize(@[0.01])
  result = params[0]
  opt.destroy()
#echo vacuumMassGivenProb(1.5e15)

proc pressureGivenEffPhotonMass(m_gamma: float, T = 4.2): float =
  ## calculates the inverse of `babyIaxoEffPhotonMass`, i.e. the pressure
  ## from a given effective photon mass in BabyIAXO
  result = m_gamma * m_gamma * T / 0.01988

template fwhmTail(name: untyped, op: untyped): untyped =
  proc `name`(m_a: float, outputInfo = false,
              massRange = none[tuple[low, high: float]]()):
    tuple[m_as, probs: seq[float], mAtFwhm, pHalf: float] =
    # calculate pressure corresponding to `m_a`
    let pressure = pressureGivenEffPhotonMass(m_a)
    # get gamma and m_gamma from pressure
    let (gamma, m_gamma) = pressureToGammaMGamma(pressure)
    doAssert abs(m_a - m_gamma) < 1e-4
    let (m_as, probs) = calcConvProbCurve(gamma, m_gamma,
                                          pressure, massRange = massRange)
    let (p, pResI) = fitToConvProb(gamma, m_gamma, pressure,
                                   nameSuffix = "teststest", createPlot = true)
    let fwhmVal = fwhm(pResI[2])
    let mAtFwhm = `op`(pResI[1], fwhmVal / 2.0)
    let pHalf = pResI[0] / 2.0
    result = (m_as: m_as, probs: probs, mAtFwhm: mAtFwhm, pHalf: pHalf)
    if outputInfo:
      echo "Axion mass m_a = ", m_a, " eV"
      echo "Pressure P = ", pressure, " mbar"
      echo "Attenuation length Gamma = ", gamma
      echo "Effective photon mass m_gamma = ", m_gamma, " eV"
      echo "FWHM ", fwhmVal

fwhmTail(lhsFwhmTail, `-`)
fwhmTail(rhsFwhmTail, `+`)

proc massDifference(m_a, prob: float, cmpMass = none[float]()): float =
  ## Given a conversion probability value at a certain mass, calculate the
  ## mass difference between this value and the vacuum curve
  var ma_cmp: float
  if cmpMass.isNone:
    ma_cmp = vacuumMassGivenProb(prob)
  else:
    let (m_as, probs, mAtFwhm, pHalf) = rhsFwhmTail(cmpMass.get)
    ma_cmp = mAtFwhm
    echo "using ", cmpMass
  echo "Corresponding ma in vacuum ", ma_cmp, " for ", prob
  result = abs(m_a - ma_cmp)

#proc gasStepsPlot(massLhs, massRhs: tuple[name: string, m: float], title = "",
#                  filterMass = none[tuple[low, high: float]]()) =
#  ## generates a plot showing the resonance curves for the two masses
#  ## (or one if one mass is 0.0). Compares it to the vacuum conversion
#  ## probability.
#  let (m_as, probs, mAtFwhm, pHalf) = lhsFwhmTail(massLhs.m, massRange = filterMass)
#  let dfLhs = seqsToDf({ "m_a / eV" : m_as,
#                         "P_a->gamma" : probs,
#  })
#  echo "LHS = ", massLhs
#  echo "LHS P = ", pressureGivenEffPhotonMass(massLhs.m)
#  var dfComb: DataFrame
#  if massRhs.name.len == 0:
#    dfComb = bind_rows([("vacuum", getDfVac(filterMass)),
#                        (massLhs.name, dfLhs)],
#                        #("xline", dfDummy),
#                        #("yline", dfDummy2),
#                       id = "case")
#  else:
#    echo "RHS = ", massRhs
#    echo "RHS P = ", pressureGivenEffPhotonMass(massRhs.m)
#    let (m_Rhs, probsRhs, _, _) = rhsFwhmTail(massRhs.m, massRange = filterMass)
#    let dfRhs = seqsToDf({ "m_a / eV" : m_Rhs,
#                           "P_a->gamma" : probsRhs,
#    })
#    dfComb = bind_rows([("vacuum", getDfVac(filterMass)),
#                        (massRhs.name, dfRhs),
#                        (massLhs.name, dfLhs)],
#                       id = "case")
#  if filterMass.isSome:
#    let fm = filterMass.get
#    dfComb = dfComb.filter(f{"m_a / eV" >= fm.low and "m_a / eV" <= fm.high})
#  var plt = ggplot(dfComb, aes("m_a / eV", "P_a->gamma", color = "case")) +
#    geom_line() +
#    scale_y_log10() +
#    ggtitle(title)
#  if filterMass.isNone:
#    plt = plt + scale_x_log10()
#  let n1 = massLhs.name.replace(" ", "_")
#  let n2 = massRhs.name.replace(" ", "_")
#  plt + ggsave(&"vacuum_helium_cutoff_{mAtFWHM}_{n1}_{n2}.pdf")


#proc massDifferenceAtMass(m_a: float, createPlot = true,
#                          outputInfo = false, cmpMass = none[float]()): float =
#  ## returns the pressure that corresponds to the mass that is required to
#  ## achieve a FWHM overlap between the vacuum curve and the buffer gas setup
#  # use `m_a` as start
#  # calculate buffer curve with this value
#  echo "Eff photon mass ", m_a
#  let (m_as, probs, mAtFwhm, pHalf) = lhsFwhmTail(m_a, outputInfo)
#  # find m_a of `pHalf` on vacuum curve
#  result = massDifference(mAtFwhm, pHalf, cmpMass)
#  if createPlot:
#    if cmpMass.isNone:
#      gasStepsPlot((name: "helium", m: m_a), massRhs = (name: "", m: 0.0))
#    else:
#      gasStepsPlot((name: "helium", m: m_a),
#                   (name: "heliumRef", m: cmpMass.get))
#
#proc findMassAtFwhm(m_a: float, cmpWithMassless = true): float =
#  ## performs minimization of `massDifferenceAtMass` to find the mass `m_a` at which
#  ## the FWHM of the buffer phase matches the vacuum line
#  proc findMin(p0: seq[float], dummy: MassFit): float =
#    if cmpWithMassless:
#      result = massDifferenceAtMass(p0[0])
#    else:
#      result = massDifferenceAtMass(p0[0], cmpMass = some(m_a))
#
#  var opt = newNloptOpt[MassFit](LN_COBYLA, 1, bounds = @[(1e-4, Inf)])
#  # hand the function to fit as well as the data object we need in it
#  var varStr = newVarStruct(findMin, MassFit())
#  opt.setFunction(varStr)
#  opt.xtol_rel = 1e-8
#  opt.ftol_rel = 1e-8
#  let (params, minVal) = opt.optimize(@[m_a])
#  result = params[0]
#  opt.destroy()
#
#when false:
#  let mstep1 = findMassAtFwhm(babyIaxoVacuumMassLimit)
#  let mstep2 = findMassAtFwhm(mstep1, cmpWithMassless = false)
#
#when false:
#  echo "\n\n------------------------------\n\n"
#  proc toPressureDiff(m1, m2: float, T = 4.2): float
#  let pdiffVacStep1 = toPressureDiff(mstep1, babyIaxoVacuumMassLimit)
#  gasStepsPlot((name: "1st He step", m: mstep1),
#               massRhs = (name: "", m: 0.0),
#               title = &"First buffer gas step; match at FWHM, ΔP = {pdiffVacStep1:.4f} mbar (@ 4.2 K)")
#  let pdiffFirst2Steps = toPressureDiff(mstep1, mstep2)
#  gasStepsPlot((name: "2nd He step", m: mstep2),
#               (name: "1st He step", m: mstep1),
#               title = &"First two He steps, each match at FWHM, ΔP = {pdiffFirst2Steps:.4f} mbar (@ 4.2 K)")
#
#proc toPressureDiff(m1, m2: float, T = 4.2): float =
#  ## assuming two masses `m1` and `m2` in `eV`, converts the
#  ## mass difference to a corresponding pressure difference
#  ## at temperature `T`
#  result = abs(pressureGivenEffPhotonMass(m2, T = T) -
#               pressureGivenEffPhotonMass(m1, T = T))
#
#
## 1 bar
#when false:
#  let mstep1Bar = findMassAtFwhm(m_gamma_1bar, cmpWithMassless = false)
#  let pdiff1bar = toPressureDiff(m_gamma_1bar, mstep1Bar)
#  gasStepsPlot((name: "Step 1", m: m_gamma_1_bar),
#               (name: "Step 2", m: mstep1Bar),
#               title = &"Two steps at 1 bar (@ 293 K) pressure, ΔP = {pdiff1bar:.4f} mbar (@ 4.2 K)",
#               filterMass = some((low: 0.22, high: 0.2999)))
#  #
#  ### 3 bar
#  let mstep3Bar = findMassAtFwhm(m_gamma_3bar, cmpWithMassless = false)
#  let pdiff3bar = toPressureDiff(m_gamma_3bar, mstep3Bar)
#  gasStepsPlot((name: "Step 1", m: m_gamma_3_bar),
#               (name: "Step 2", m: mstep3Bar),
#               title = &"Two steps at 3 bar (@ 293 K) pressure, ΔP = {pdiff3bar:.4f} mbar (@ 4.2 K)",
#               filterMass = some((low: 0.42, high: 0.48)))
#
#proc m_a_vs_density(pstart, pstop: float) =
#  let pressures = logspace(pstart.log10, pstop.log10, 1000)
#  let densities = pressures.mapIt(density(it, 4.2))
#  let masses = pressures.mapIt(babyIaxoEffMass(it)) # corresponding masses
#  # convert both seqs to a dataframe
#  let ref1Bar = density(1000, 293.15)
#  let df1bar = seqsToDf({"ρ / g/cm^3" : @[ref1Bar, ref1Bar], "m_a / eV" : @[1e-2, 1.0]})
#  let ref3Bar = density(3000, 293.15)
#  let df3bar = seqsToDf({"ρ / g/cm^3" : @[ref3Bar, ref3Bar], "m_a / eV" : @[1e-2, 1.0]})
#  let refVacLimit = density(pressureGivenEffPhotonMass(babyIaxoVacuumMassLimit))
#  let dfVacLimit = seqsToDf({"ρ / g/cm^3" : @[refVacLimit, refVacLimit], "m_a / eV" : @[1e-2, 1.0]})
#  let df = seqsToDf({"ρ / g/cm^3" : densities, "m_a / eV" : masses})
#  let dfComb = bind_rows([("ma vs ρ", df),
#                          ("1 bar @ 293 K", df1bar),
#                          ("3 bar @ 293 K", df3bar),
#                          ("Vacuum limit", dfVacLimit)],
#                         id = "Legend")
#  ggplot(dfComb, aes("ρ / g/cm^3", "m_a / eV", color = "Legend")) +
#    geom_line() +
#    scale_x_log10() +
#    scale_y_log10() +
#    ggtitle("Sensitive axion mass in eV depending on helium density in g / cm^3") +
#    ggsave("m_a_vs_density.pdf")

when false:
  m_a_vs_density(pressureGivenEffPhotonMass(babyIaxoVacuumMassLimit) * 0.9,
                 pressureGivenEffPhotonMass(mstep3Bar) * 1.1)
