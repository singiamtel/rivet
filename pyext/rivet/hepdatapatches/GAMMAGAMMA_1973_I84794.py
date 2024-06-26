import yoda

def patch(path, ao):
    # reorder points
    if "/REF/GAMMAGAMMA_1973_I84794" in path and ("d03" in path or "d04" in path) :
        edges = ao.xEdges()
        edges.reverse()
        newao = yoda.BinnedEstimate1D(edges, ao.path())
        for val in edges :
            b = ao.binAt(val)
            newao.set(newao.binAt(val).index(),b)
        ao = newao
    return ao
