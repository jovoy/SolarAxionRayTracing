import std / [math, strutils, algorithm, random, sequtils, os, strformat, tables, sugar, strscans]
import unchained

proc rij_s(ni: float, nj: float, theta0: Degree, n0: float): float =
  let a = n0 * cos(theta0.to(Radian))
  let x = ni * sqrt(1.0 - pow((a/ni), 2.0))
  let y = nj * sqrt(1.0 - pow((a/nj), 2.0))
  result = ((x-y)/(x+y))

proc rij_p(ni: float, nj: float, theta0: Degree, n0: float): float =
  let a = n0 * cos(theta0.to(Radian))
  let x = nj * sqrt(1.0 - pow((a/ni), 2.0))
  let y = ni * sqrt(1.0 - pow((a/nj), 2.0))
  result = ((x-y)/(x+y))

proc fi(lambd: float, di: float, ni: float, theta0: Degree, n0: float): float =
  result = (4.0 * PI /lambd) * di * ni *  sqrt(1.0 - pow((n0 * cos(theta0.to(Radian)) /ni), 2.0))