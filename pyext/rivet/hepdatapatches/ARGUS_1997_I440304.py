import yoda
def patch(path, ao):
    if "ARGUS_1997_I440304" in path and "d01" in path:
        newHist1 = yoda.BinnedEstimate1D([0.7,1.0],path)
        newHist1.setAnnotation("IsRef",1)
        newHist1.set(1,ao.bin(2))
        newHist2 = yoda.BinnedEstimate1D([0.0,1.0],path.replace("x01","x02"))
        newHist2.setAnnotation("IsRef",1)
        newHist2.set(1,ao.bin(1))
        return [newHist1,newHist2]
    return ao
