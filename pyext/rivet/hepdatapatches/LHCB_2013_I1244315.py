import yoda
def patch(path, ao):
    # fix hist, this one is a mess really 3 different 2d histos
    if ("LHCB_2013_I1244315" in path and ("d01" in path or "d02" in path)) :
        edges=ao.xEdges()
        newHists=[]
        newHist = yoda.BinnedEstimate1D(edges,ao.path())
        for i in range(1,4) :
            path = ao.path().replace("x01","x0%s"%i)
            newHists.append(yoda.BinnedEstimate1D(edges,path))
            newHists[-1].setTitle(ao.title())
            newHists[-1].setAnnotation("IsRef",1)
        for iy in range(1, ao.numBinsY()+1):
            for ix in range(1, ao.numBinsX()+1):
                b = ao.bin(ix,iy)
                newHists[iy-1].set(ix, b)
        return newHists
    else :
        return ao

