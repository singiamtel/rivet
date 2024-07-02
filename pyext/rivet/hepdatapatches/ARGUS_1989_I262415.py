import yoda
def patch(path, ao):
    if "ARGUS_1989_I262415" in path and ("d01" in path or "d02" in path):
        newHist1=yoda.BinnedEstimate1D([0.1,0.4],path)
        newHist1.set(1,ao.bin(2))
        newHist2=yoda.BinnedEstimate1D([0.,1.],path.replace("x01","x02"))
        newHist2.set(1,ao.bin(3))
        return [newHist1,newHist2]
    return ao
