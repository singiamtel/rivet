import yoda

def patch(path, ao):
    if "CMS_2017_I1631985/d03-x01-y01" in path or "CMS_2017_I1631985/d04-x01-y01" in path:
        newaos = [ ao ]
        h = ao.path()
        hbase = '/'.join(h.split('/')[:-1])
        hbase += '/d09-x01-y0' if 'd03' in h else '/d10-x01-y0'
        for j in range(6):
            hname = hbase + '%i' % (j+1)
            newao = yoda.Estimate1D(ao.xEdges(), hname)
            for i in range(1, newao.numBins()+1):
                b = ao.bin(i,j+1)
                newao.set(i, b)
            newaos.append(newao)
        return newaos
    return ao
