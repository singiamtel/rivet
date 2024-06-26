import yoda
def patch(path, ao):
    if "DELPHI_2000_I524693" in path and "d01" in path:
        edges = ao.xEdges()
        edges.append(206)
        newAO = yoda.BinnedEstimate1D(edges,path)
        for i in range(0,ao.numBins()) :
            newAO.set(i+1,ao.bin(i+1))
        if "y01" in path :
            newAO.bin(ao.numBins()+1).set(28.72,0.70)
        elif "y02" in path :
            newAO.bin(ao.numBins()+1).set(36.32,3.08)
        elif "y03" in path :
            newAO.bin(ao.numBins()+1).set(24.22,1.09)
        elif "y04" in path :
            newAO.bin(ao.numBins()+1).set(4.50,1.05)
        return newAO
    return ao


