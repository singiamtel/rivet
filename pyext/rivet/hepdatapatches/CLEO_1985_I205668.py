import yoda
def patch(path, ao):
    if "CLEO_1985_I205668" in path and "d01-x01-y02" in path:
        # remove extra point
        edges = ao.xEdges()
        edges.remove("0.105")
        print(edges)
        aonew=yoda.BinnedEstimate1D(edges,ao.path())
        ioff=0
        for i in range(1,len(edges)+1) :
            if i==4 : ioff+=1
            aonew.set(i,ao.bin(i+ioff))
        return aonew
    return ao
