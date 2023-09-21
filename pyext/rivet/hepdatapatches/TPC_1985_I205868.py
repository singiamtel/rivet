import yoda
def patch(path, ao):
    # fix type of hist from 2D to 1D, 2 equivalent independent axes
    if "TPC_1985_I205868" in path:
        newao = yoda.BinnedEstimate1D(ao.xEdges(), ao.path())
        for i in range(1, newao.numBins()+1):
            b = ao.bin(i,i)
            newao.set(i, b)
        ao = newao
    return ao
