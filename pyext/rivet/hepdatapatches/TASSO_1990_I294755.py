import yoda
def patch(path, ao):
    # remove average bins
    if "TASSO_1990_I294755" in path and "d05" in path and ("y01" in path or "y02" in path ):
        newao = yoda.BinnedEstimate1D(ao.xEdges()[:-1],path)
        for i in range(0,ao.numBins()-1) :
            newao.set(i+1,ao.bin(i+1))
        return newao
    return ao
