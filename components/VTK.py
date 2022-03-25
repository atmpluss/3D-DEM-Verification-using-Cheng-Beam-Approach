def write_vtk_file(particleData,step, folder, factor=1):
    outputStr = "# vtk DataFile Version 3.0\nParticle data\nASCII\nDATASET POLYDATA\n\n"
    numParticles=len(particleData)
    outputStr += "POINTS " + str(numParticles) + " float\n"
    for p in particleData:
        outputStr += str(p.pos[0]) + " " + str(p.pos[1]* factor) + " " + str(p.pos[2]) + "\n"
    outputStr += "\nPOINT_DATA " + str(numParticles) + "\n"
    outputStr += "SCALARS radius float\nLOOKUP_TABLE default\n"
    for p in particleData:
        outputStr += str(p.r) + "\n"
    outputStr +="\nVECTORS velocity float\n"
    for p in particleData:
        outputStr += str(p.v[0]) + " " + str(p.v[1]) + " " + str(p.v[2]) + "\n"
    with open(str(folder)+"/particles"+str(step)+".vtk","w") as fout:
        fout.write(outputStr)
