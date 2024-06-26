import yoda
def patch(path, ao):
    # fix bin widths
    if ("SLD_1999_I469925" in path and
        ("d01" in path or "d02" in path or "d03" in path or
         "d05" in path or "d06" in path or "d07" in path or
         "d08" in path or "d09" in path)) :
        edges=ao.xEdges()
        edge1=[]
        edge2=[]
        masked=[]
        emax = len(edges)-1
        if ("d01-x01-y01" in path or "d02-x01-y01" in path or
            "d03-x01-y01" in path)  : emax+=1
        for i in range(0,emax) :
            tmp=edges[i].split("-")
            tmp[0]=tmp[0].strip()
            tmp[1]=tmp[1].strip()
            if len(edge1) == 0 :
                edge1.append(tmp[0])
            elif tmp[0] !=edge1[-1] :
                masked.append(len(edge1))
                edge1.append(tmp[0])
            edge1.append(tmp[1])
        for i in range(0,len(edge1)) : edge1[i]=float(edge1[i])
        ao1=yoda.BinnedEstimate1D(edge1,path)
        if len(masked)!=0 : ao1.maskBins(masked)
        j=1
        for i in range(1, ao1.numBins()+1):
            b = ao.bin(i)
            if i in masked : j+=1
            ao1.set(j, b)
            j+=1
        if "d01" not in path and "d02" not in path and "d03" not in path:
            tmp=edges[-1].split("-")
            tmp[0]=tmp[0].strip()
            tmp[1]=tmp[1].strip()
            edge2=[float(tmp[0]),float(tmp[1])]
            ao2=yoda.BinnedEstimate1D(edge2,path.replace("x01","x02"))
            ao2.set(1,ao.bin(ao.numBins()))
            return [ao1,ao2]
        else:
            return ao1
    return ao
