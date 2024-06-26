import yoda

def patch(path, ao):
    # fix hist, really a 2D histo with two equivalent x axes not a 3d one
    if "JADE_1990_I282847" in path and ("d03" in path or "d04" in path or "d05" in path):
        sqrtS= 44. if 'd04' in path else 35.
        newedges = [ 2.* x / sqrtS for x in ao.xEdges() ]
        newao = yoda.BinnedEstimate1D(newedges, ao.path())
        for i in range(1, newao.numBins()+1):
            b = ao.bin(i,i)
            newao.set(i, b)
        ao = newao
    return ao
