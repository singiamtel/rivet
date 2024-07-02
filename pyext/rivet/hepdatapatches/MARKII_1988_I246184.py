
import yoda
def patch(path, ao):
    # last bin not properly normalised (strip it off)
    if "MARKII_1988_I246184" in path and ("d11" in path or "d12" in path or "d14" in path or
                                          "d29" in path or "d30" in path or "d32" in path or
                                          "d47" in path or "d48" in path or "d50" in path ) :
        newao=yoda.BinnedEstimate1D(ao.xEdges()[:-1],ao.path())
        for i in range(0,newao.numBins()) :
            newao.set(i+1,ao.bin(i+1))
        return newao
    return ao

