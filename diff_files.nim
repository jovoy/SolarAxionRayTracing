import os

let fIg = "/home/oy/Documents/GitHub/AxionElectronLimit/OPCD_3.3/mono/fm01.mesh"
let f1 = readFile("/home/oy/Documents/GitHub/AxionElectronLimit/OPCD_3.3/mono/fm02.mesh")
for file in walkFiles("/home/oy/Documents/GitHub/AxionElectronLimit/OPCD_3.3/mono/fm??.mesh"):
  if file != fIg:
    echo "Comparing file ", file
    let f2 = readFile(file)
    doAssert f1 == f2, "file differs from " & $file
