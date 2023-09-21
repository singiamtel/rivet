import yoda

def patch(path, ao):
    if path == '/REF/ATLAS_2018_I1711223/d12-x01-y01':
        edges = [ '2', '3', '4', r'$\geq5$' ]
        newao = yoda.BinnedEstimate1D(edges, path)
        for b in ao.bins():
            newao.set(b.index(), b)
        ao = newao
    if path == '/REF/ATLAS_2018_I1711223/d20-x01-y01':
        edges = [ '0', '1', '2', r'$\geq3$' ]
        newao = yoda.BinnedEstimate1D(edges, path)
        for b in ao.bins():
            newao.set(b.index(), b)
        ao = newao
    return ao

