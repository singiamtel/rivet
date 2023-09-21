import yoda, math

def patch(path, ao):
    if "HRS_1987_I215848" in path and ("d02" in path or "d03" in path or "d04" in path):
        # bins edges given as <z>, convert momentum limits to energy fractions instead
        # momentum limits specified in paper:
        edges = [0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9]
        ix=0
        mass=0.
        if "d02" in path:
            mass =0.13957
        elif "d03" in path:
            mass=0.493677
        elif "d04" in path:
            mass=0.938272
        edges = [ 2.*math.sqrt(e**2+mass**2)/29. for e in edges ]
        newao = yoda.BinnedEstimate1D(edges, ao.path())
        for b in ao.bins():
            newao.set(b.index(), b)
    return ao
