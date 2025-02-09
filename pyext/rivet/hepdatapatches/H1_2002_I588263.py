import math

def patch(path, ao):
    # remove dummy axes
    if 'd01' in path or 'd02' in path:  return ao
    newao = yoda.BinnedEstimate1D(ao.xEdges(), path)
    idx = 1
    for b in ao.bins():
        if math.isnan(b.val()):  continue
        newao.set(idx, b)
        idx += 1
    return newao

