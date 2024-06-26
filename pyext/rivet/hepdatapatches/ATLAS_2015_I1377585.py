import yoda,math

def patch(path, ao):
    if "ATLAS_2015_I1377585" in path and "d01" in path:
        # 3D hist is really only 1D hist
        newAO = yoda.BinnedEstimate1D([7000],path)
        for b in ao.bins() :
            if not math.isnan(b.val()) :
                newAO.set(1,b)
        return newAO
    return ao
