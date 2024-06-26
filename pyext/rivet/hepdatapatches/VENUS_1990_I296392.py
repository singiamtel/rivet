import yoda

def patch(path, ao):
    # fix bin widths (the ones in hepdata are for p not x)
    if path == "/REF/VENUS_1990_I296392/d02-x01-y01" :
        bins =  [ 2.*b/34.4 for b in ao.xEdges() ]
        newao = yoda.BinnedEstimate1D(bins, ao.path())
        for b in ao.bins():
            newao.set(b.index(), b)
        ao = newao
    return ao
