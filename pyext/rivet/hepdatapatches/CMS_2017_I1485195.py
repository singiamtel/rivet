import yoda,math
def patch(path, ao):
    # make 3d scatters 2d
    if "CMS_2017_I1485195" in path and type(ao) == yoda.core.BinnedEstimate2D :
        if "d05" in path : 
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
        elif "d06" in path :
            edges=ao.yEdges()
            newHist = yoda.Estimate1D(edges,ao.path())
            newHist.setTitle(ao.title())
            newHist.setAnnotation("IsRef",1)
            for b in ao.bins() :
                if math.isnan(b.val()) : continue
                for i in range(1,len(edges)) :
                    if b.yMid() >edges[i-1] and b.yMid() < edges[i] :
                        newHist.set(i,b)
            return newHist
    return ao
