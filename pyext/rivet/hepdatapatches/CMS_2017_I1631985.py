
def patch(path, ao):
    # fix bin widths
    if "CMS_2017_I1631985/d03-x01-y01" in path or "CMS_2017_I1631985/d04-x01-y01" in path:
        newaos = [ ao ]
        hbase = '/'.join(h.split('/')[:-1])
        hbase += '/d09-x01-y0' if 'd03' in h else '/d10-x01-y0'
        for j in range(6):
            hname = hbase + '%i' % (j+1)
            newao = yoda.Scatter2D(hname)
            for i in range(34):
                p = ao.point(6*i + j)
                newao.addPoint(p.x(), p.z(), p.xErrs(), p.zErrs())
            newaos.append(newao)
        return newaos
    return ao
