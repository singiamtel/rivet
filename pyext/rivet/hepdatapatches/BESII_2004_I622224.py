import yoda,math

# removal of average bins and replace bin no with x values
def patch(path, ao):
    if "BESII_2004_I622224" in path :
        ioff=0
        if   "d10" in path : ioff=1
        elif "d09" in path : ioff=2
        elif "d08" in path : ioff=3
        else : return ao
        edges = ao.xEdges()[ioff:]
        newao = yoda.BinnedEstimate1D(edges, ao.path())
        for i in range(ioff,ao.numBins()) :
            newao.set(i+1-ioff,ao.bin(i+1))
        return newao
    return ao
