import std / [sequtils, strformat, algorithm, options, strutils]
import seqmath, ggplotnim, nlopt, unchained

proc density(p: float, temp:float ): float =
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

proc numDensity(c: float): float =
  ## converts a molar concentration in mol / m^3 to a number density
  const Z = 2 # number of electron in helium atom
  result = Z * 6.022e23 * c

proc effPhotonMass*(ne: float): float =
  ## returns the effective photon mass for a given electron number density
  const alpha = 1.0 / 137.0
  const me = 511e3 # 511 keV
  # note the 1.97e-7 cubed to account for the length scale in `ne`
  result = sqrt( pow(1.97e-7, 3) * 4 * PI * alpha * ne / me )

proc molarAmount(p, vol, temp: float): float =
  ## calculates the molar amount of gas given a certain pressure,
  ## volume and temperature
  ## the pressure is assumed in mbar
  const gasConstant = 8.314 # joule K^-1 mol^-1
  let pressure = p * 1e2 # pressure in Pa
  result = pressure * vol / (gasConstant * temp)
  #echo "Molar amount for P = ", pressure, " Pa is ", result

proc effPhotonMass2*(p : float, length : float, radBore : float, temp : float): float =
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

proc momentumTransfer(m_gamma, m_a: float, E_a = 4.2.keV): float =
  ## calculates the momentum transfer for a given effective photon
  ## mass `m_gamma` and axion mass `m_a` at an axion energy of
  ## 4.2 keV `E_a` (by default).
  #const c = 299792458
  result = abs((m_gamma * m_gamma - m_a * m_a) / (2 * E_a.to(eV).float))

proc logMassAttenuation(e: keV): float =
  ## calculates the logarithm of the mass attenuation coefficient for a given
  ## energy `e` in `keV` and the result in `cm^2/g`
  result = -1.5832 + 5.9195 * exp(-0.353808 * e.float) + 4.03598 * exp(-0.970557 * e.float)

proc axionConversionProb2*(m_a: float, energyAx: keV, pressure: float, temp: float, length: float, radBore: float, g_agamma: float, B: float): float =
  ## given an axion mass and an inverse absorption length
  ## calculates the conversion probabilty
  # both `g_agamma` and `B` only scale the absolute value `P`, does not matter
  let energyGamma = energyAx #in keV
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

proc intensitySuppression2*(energy: keV, distanceMagnet: float, distancePipe: float, pressure: float, tempMagnet: float, tempPipe: float64): float =
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

when isMainModule:
  echo effPhotonMass2(36.61, 10.0, 0.35, 100.0) #for 1 bar at room temperature: T=4.2K -> p=15.38mbar  and T=100K -> p=366.1mbar -> axion_mass=0.2698eV
                                                #for 0.5 bar at room temperature: T=4.2K -> p=xmbar  and T=100K -> p=183.05mbar -> axion_mass=0.1907eV
                                                #for 0.3 bar at room temperature: T=4.2K -> p=xmbar  and T=100K -> p=109.8mbar -> axion_mass=0.1477eV
                                                #for 0.1 bar at room temperature: T=4.2K -> p=xmbar  and T=100K -> p=36.61mbar -> axion_mass=0.0853eV
