import yoda,math

def patch(path, ao):
    if "BELLE_2023_I2624324/d03-x01-y01" in path or "BELLE_2023_I2624324/d06-x01-y01" in path:
        bins = [[1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.51],
                [-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00],
                [-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00],
                [0.]]
        for i in range(1,11) : bins[3].append(0.2*math.pi*i)
        if "d03" in path :
            ih=2
            jmax=2
        else :
            ih=1
            jmax=5
        newaos = [ao]
        for j in range(1,jmax) :
            for i in range(1,5) :
                newpath = "/REF/BELLE_2023_I2624324/d0%s-x0%s-y0%s"%(ih,i,j)
                newao = yoda.Estimate1D(bins[i-1], newpath)
                newao.setPath(newpath)
                newao.setAnnotation("IsRef",1)
                newaos.append(newao)
                for k in range(1,11) :
                    width = (bins[i-1][k]-bins[i-1][k-1])
                    oldb = ao.bin(40*(j-1)+10*(i-1)+k)
                    newao.bin(k).setVal(oldb.val()/width)
                    sources = oldb.sources()
                    errs = [ (s, oldb.err(s)) for s in sources ]
                    for label, epair in errs:
                        newao.bin(k).setErr(epair, label)

                    
        return newaos
    return ao
