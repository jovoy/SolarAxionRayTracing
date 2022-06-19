import std / [httpclient, strutils, os, parseutils, strformat, math, strscans, sugar, random]
import seqmath

import cligen/[osUt,mslice,procpool]

const henkeUrl = "https://henke.lbl.gov/"
const materialInterface = "cgi-bin/laymir.pl"

# Constants for the output of the downloaded files
const outpath = "../resources/henke_download/"
const fileTmpl = "$#degGold$#microns.csv"

type
  ScanKind = enum
    skEnergy = "Energy"               # from 30 eV to 30,000 eV
    skWavelength = "Wavelength"       # from 0.041 nm to 41 nm
    skAngle = "Angle"                 # from 0 to 90°
  FixedKind = enum
    fkEnergy = "Energy+%28eV%29"
    fkAngle = "Angle+%28deg%29"
  PlotKind = enum
    pkLinear = "Linear"
    pkLog = "Log"
    pkLinLog = "LinLog"
    pkLogLin = "LogLin"
  OutputKind = enum
    okPlot = "Plot"
    okEps = "EPS"
  PolarizationKind = enum
    polP = -1
    polUnpolarized = 0
    polS = 1
  HenkeReq = object
    layer: string                        # chemical formula of layer
    layerDensity: float                  # layer density in g/cm³ (if -1: tabulated value density)
    layerThickness: float                # thickness of layer in nm
    surfaceRoughness: float              # top surface roughness in nm ("sigma1")
    substrateMaterial: string            # chemical formula of the substrate material
    substrateDensity: float              # density of the substrate material (if -1: tab. value density)
    substrateRoughness: float            # top surface roughness in nm of substrate ("sigma2")
    polarization: PolarizationKind       # polarization (s polarization: 1, p polarization: -1, unpol.: 0)
    scanQuantity: ScanKind               # quantity to scan over
    scanMin: range[0.0 .. 30000.0]       # minimum value from which to scan
    scanMax: range[0.0 .. 30000.0]       # minimum value from which to scan
    numberOfPoints: range[0 .. int.high] # number of scan points (max 500; if more we split the download into multiple parts)
    fixedQuantity: FixedKind             # quantity to keep fixed
    fixedValue: range[0.0 .. 30000.0]    # fixed value of the quantity
    plotKind: PlotKind                   # kind of plot to produce
    outputKind: OutputKind

proc toStr(f: float): string =
  if f == -1.0:
    result = "-1"
  else:
    result = $f

proc `$`(h: HenkeReq): string =
  result.add "Layer=" & $h.layer & "&"
  result.add "Ldensity=" & h.layerDensity.toStr & "&"
  result.add "Thick=" & h.layerThickness.toStr & "&"
  result.add "Sigma1=" & h.surfaceRoughness.toStr & "&"
  result.add "Substrate=" & $h.substrateMaterial & "&"
  result.add "Sdensity=" & h.substrateDensity.toStr & "&"
  result.add "Sigma2=" & h.substrateRoughness.toStr & "&"
  result.add "Pol=" & $(h.polarization.ord) & "&"
  result.add "Scan=" & $h.scanQuantity & "&"
  result.add "Min=" & $h.scanMin & "&"
  result.add "Max=" & $h.scanMax & "&"
  result.add "Npts=" & $(h.numberOfPoints - 1) & "&" # henke takes number of *steps* and not data points
  result.add "temp=" & $h.fixedQuantity & "&"
  result.add "Fixed=" & $h.fixedValue & "&"
  result.add "Plot=" & $h.plotKind & "&"
  result.add "Output=" & $h.outputKind

proc initHenkeReq(layer = "Au",
                   layerDensity = -1.0,
                   layerThickness = 250.0,
                   surfaceRoughness = 0.0,
                   substrateMaterial = "SiO2",
                   substrateDensity = -1.0,
                   substrateRoughness = 0.0,
                   polarization = polUnpolarized,
                   scanQuantity = skEnergy,
                   scanMin = 30.0,
                   scanMax = 5000.0,
                   numberOfPoints = 500,
                   fixedQuantity = fkAngle,
                   fixedValue = 1.0,
                   plotKind = pkLinear,
                   outputKind = okPlot): HenkeReq =
  result = HenkeReq(layer: layer,
                     layerDensity: layerDensity,
                     layerThickness: layerThickness,
                     surfaceRoughness: surfaceRoughness,
                     substrateMaterial: substrateMaterial,
                     substrateDensity: substrateDensity,
                     substrateRoughness: substrateRoughness,
                     polarization: polarization,
                     scanQuantity: scanQuantity,
                     scanMin: scanMin,
                     scanMax: scanMax,
                     numberOfPoints: numberOfPoints,
                     fixedQuantity: fixedQuantity,
                     fixedValue: fixedValue,
                     plotKind: plotKind,
                     outputKind: outputKind)

proc testHenkeStr() =
  let testExp = "Layer=Au&Ldensity=-1&Thick=250.0&Sigma1=0.0&Substrate=SiO2&Sdensity=-1&Sigma2=0.0&Pol=0&Scan=Energy&Min=7495.0&Max=15000.0&Npts=499&temp=Angle+%28deg%29&Fixed=0.95&Plot=LinLog&Output=Plot"
  let req = initHenkeReq(fixedValue = 0.95,
                          scanMin = 7495.0,
                          scanMax = 15000.0,
                          plotKind = pkLinLog)
  doAssert testExp == $req
testHenkeStr()

proc postHenke(client: HttpClient, h: HenkeReq): string =
  if h.numberOfPoints > 500:
    raise newException(ValueError, "Number of points to download may not exceed 500!")
  let header = newHttpHeaders({"Content-Type" : "application/x-www-form-urlencoded"})
  let postPath = henkeUrl & materialInterface
  let resp = client.request(postPath, httpMethod = HttpPost, body = $h)
  if resp.status == $Http200:
    result = resp.body
  else:
    raise newException(HttpRequestError, "Request to henke with " & $h & " failed!")

proc extractDataPath(resp: string): string =
  ## Extrats the path from the response body that contains the requested data file
  for l in resp.splitLines:
    let findIt = "<a HREF=\""
    let idx = find(l, findIt)
    if idx > 0:
      let num = l.parseUntil(result, '\"', findIt.len + idx)
      break

proc downloadDatafile(client: HttpClient, file: string, angle: float): seq[string] =
  ## Downloads the data file from the Henke server and returns is split by lines
  ## for easier handling of the header.
  echo "Downloading the path: ", file, " for angle ", angle
  var data = client.getContent(henkeUrl & file)
  result = data.strip.splitLines
  var idx = 0
  for i, l in mpairs(result):
    if i > 1: break
    if i == 0: l = "# " & l.strip
    else: l = "#" & l.strip.replace(" ", "").replace(",", " ")
  var l0 = result[0]
  result[0] = result[1]
  result[1] = l0

proc storeDatafile(data: seq[string], h: HenkeReq) =
  let angle = &"{h.fixedValue.float:.6f}"
  let thickness = &"{h.layerThickness / 1000.0:.2f}"
  let outfile = outpath / fileTmpl % [angle, thickness]
  createDir(outpath)
  let dataStr = data.join("\n")
  writeFile(outfile, dataStr.strip & "\n")

proc sanityCheckHeader(h: HenkeReq, data: seq[string], path: string): bool =
  let header = data[1] # care only about line 1 similar to:
  # Au 250.nm on SiO2 at 1.3363deg, P=0.
  let prefix = "# Au 250.nm on SiO2 at "
  let (success, angle) = scanTuple(header.dup(removePrefix(prefix)), "$fdeg")
  if not success:
    result = false
    echo "ERROR: Could not parse angle from data header. Header was: " & $header &
      " and post: " & $h & " from path " & path
  else:
    result = abs(angle - h.fixedValue) < 1e-4
    if not result:
      echo "ERROR: Expected " & $h.fixedValue & " but got " & $angle &
        ". Something broke! We broke henke maybe. The line " & $header & " and post " & $h & " from path " & $path & "\n\n\n\n\n"

proc download(client: HttpClient, h: HenkeReq): bool =
  ## Takes care of downloading (all) the file(s) for the given Henke request.
  ## This means if there are more than 500 data points requested, it will be
  ## split into multiple data files
  ##
  ## Returns `true` if the file was downloaded successfully.
  const maxPoints = 500
  let scanVals = linspace(h.scanMin, h.scanMax, h.numberOfPoints) # all values to scan
  let nchunks = ceilDiv(h.numberOfPoints, maxPoints)
  var data: seq[string]
  var pointsLeft = h.numberOfPoints
  for chunk in 0 ..< nchunks:
    var hl = h
    hl.numberOfPoints = min(pointsLeft, maxPoints) # remaining points or max points
    hl.scanMin = scanVals[chunk * maxPoints]       # get min and max from values
    hl.scanMax = scanVals[min((chunk + 1) * maxPoints, h.numberOfPoints) - 1]
    echo "INFO: Downloading chunk ", chunk, " in range ", hl.scanMin, " ↦ ", hl.scanMax, " for ", hl.fixedValue
    let dataPath = client.postHenke(hl).extractDataPath()
    let chkData = client.downloadDatafile(dataPath, h.fixedValue)
    if not h.sanityCheckHeader(chkData, dataPath):
      return false # file did not pass sanity check. Abort and try again
    if data.len == 0: # simply add all data including header
      data.add chkData
    else:             # only add data without header
      data.add chkData[2 ..< chkData.len]
    if chunk < nchunks - 1: # in last chunk cannot decrease, as `pointsLeft` is a positive range
      dec pointsLeft, maxPoints
  if h.numberOfPoints != data.len - 2:
    result = false
  else:
    storeDatafile(data, h)
    result = true

proc downloadAllAngles(client: HttpClient,
                       energyMin = 30.0,
                       energyMax = 15000.0,
                       angleMin = 0.05,
                       angleMax = 1.0,
                       numAngles = 100,
                       numEnergies = 1000,
                       jobs = 8) =
  let angles = linspace(angleMin, angleMax, numAngles)
  var pp = initProcPool((
    proc(r, w: cint) =
      let i = open(r)
      var angle: float
      while i.uRd(angle):
        if jobs > 1:
          var rnd = initRand((angle * 10000.0).round.int)
          sleep(rnd.rand(jobs * 100)) # randomly sleep a bit to space out the requests to decrease likelihood
                                      # of bad server responses (more or less seems to work)
        let req = initHenkeReq(fixedValue = angle,
                               numberOfPoints = numEnergies,
                               scanMin = energyMin,
                               scanMax = energyMax)
        echo "INFO: Requesting angle: ", angle
        var success = false
        while not success:
          # failed once, try once more. If we fail again, stop
          try:
            success = client.download(req)
            if not success:
              echo "WARN: Failed to download due to missing data at: ", $angle, ". Trying again..."
          except HttpRequestError:
            echo "WARN: Failed to download due to HttpRequestError: ", $req, ". Trying again..."

        discard w.wrLine "INFO: Finished file angle: " & $angle
    ), framesLines, jobs)
  proc prn(s: MSlice) = echo s
  pp.evalOb angles, prn

proc main(jobs: int) =
  let client = newHttpClient()

  let ϑs = linspace(0.0, 1.5, 1000) # [0.5] # the angle under which we check
  client.downloadAllAngles(energyMin = 30.0, energyMax = 15000.0,
                           angleMin = 0.0, angleMax = 1.5,
                           numAngles = 1000, numEnergies = 1000,
                           jobs = jobs)

when isMainModule:
  import cligen
  dispatch main

when false:
  import datamancer
  block TestParse:
    let df = readCsv("/home/oy/Documents/GitHub/AxionElectronLimit/resources/henke_download/1.00degGold0.25microns.csv", sep = ' ', header = "#")
    echo df
    doAssert df.len == 500
    doAssert df.getKeys().len == 3
