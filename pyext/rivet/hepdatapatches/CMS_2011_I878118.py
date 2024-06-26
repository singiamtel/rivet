import yoda,math
def patch(path, ao):
    # fix hist, really a 1D histo with average of x values added
    if ("CMS_2011_I878118" in path and
        ("d05" in path or "d06" in path or "d07" in path)) :
        edges=ao.xEdges()
        newHist = yoda.BinnedEstimate1D(edges,ao.path())
        newHist.setTitle(ao.title())
        newHist.setAnnotation("IsRef",1)
        for bin in ao.bins() :
            if math.isnan(bin.val()) : continue
            for i in range(1,len(edges)) :
                if bin.xMid() >edges[i-1] and bin.xMid() < edges[i] :
                    newHist.set(i,bin)
        return newHist
    return ao
