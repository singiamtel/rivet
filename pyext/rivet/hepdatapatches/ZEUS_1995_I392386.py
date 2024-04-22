import yoda

def patch(path, ao):
    # fix type of hist from 2D to 1D, 2 equivalent independent axes
    if "ZEUS_1995_I392386" in path:
        if any([ d in path for d in [ 'd01', 'd02', 'd03', 'd04', 'd05', 'd06', 'd07', 'd08' ]):
            newao = yoda.BinnedEstimate1D(ao.xEdges(), ao.path())
            for i in range(1, newao.numBins()+1):
                b = ao.bin(i,i)
                newao.set(i, b)
            ao = newao
    return ao
