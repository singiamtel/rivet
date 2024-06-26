import yoda,math
# slice 2d objects into 1d objects
def patch(path, ao):
    if "BESIII_2023_I2633025" in path:
        tmp = '/REF/BESIII_2023_I2633025/d%02d-x01-y%02d'
        n = 3 if "d01" in path else 4
        out = [ ao ]
        bins=[-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
        total=0.
        for b in ao.bins() :
            if not math.isnan(b.val()) : total+=b.val()
        temp=[]
        for val in ao.yEdges() : temp.append(float(val))
        temp.sort()
        for i,val in enumerate(temp) :
            iy = 1+ao.yEdges().index(str(val))
            newpath = tmp % (n, i+1)
            imin=0
            imax=20
            for ix in range(1, ao.numBinsX()+1):
                if math.isnan(ao.bin(ix,iy).val()) :
                    if ix<=10 : imin=max(imin,ix)
                    else      : imax=min(imax,ix)
            newao = yoda.BinnedEstimate1D(bins[imin:imax], newpath)
            for ix in range(1, ao.numBinsX()+1):
                b = ao.bin(ix,iy)
                if math.isnan(b.val()) : continue
                newao.set(ix-imin, b)
            newao.scale(1./total/1e-2)
            out.append(newao)
        return out
    return ao
