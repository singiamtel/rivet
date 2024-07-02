import yoda,math
def patch(path, ao):
    # fix hist, really a 2D histo with average of x values added
    if ("CMS_2015_I1342266" in path and
        ("d01" in path or "d02" in path or "d03" in path)) :
        edges=ao.xEdges()
        newHist = yoda.Estimate1D(edges,ao.path())
        newHist.setTitle(ao.title())
        newHist.setAnnotation("IsRef",1)
        for b in ao.bins() :
            if math.isnan(b.val()) : continue
            for i in range(1,len(edges)) :
                if b.xMid() >edges[i-1] and b.xMid() < edges[i] :
                    newHist.set(i,b)
        return newHist
    return ao
