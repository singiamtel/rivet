import yoda

def patch(path, ao):
    # 2D hist is really only 1D hist
    if "ARGUS_1990_I295621" in path and ("d07" in path or "d08" in path):
        newAO = yoda.BinnedEstimate1D(["$\\Upsilon_{1.2}/\\textrm{continuum}$"],path)
        newAO.set(1,ao.bin(3))
        return newAO
    return ao
