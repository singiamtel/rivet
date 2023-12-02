import yoda

def patch(path, ao):
    needs_patching = [
        '/REF/ATLAS_2022_I2614196/d15-x01-y01',
        '/REF/ATLAS_2022_I2614196/d16-x01-y01',
    ]
    if path in needs_patching:
        tmp = '/REF/ATLAS_2022_I2614196/d%d-x01-y0%d'
        n = 51 if 'd15' in path else 52
        out = [ ao ]
        for iy in range(1, ao.numBinsY()+1):
            newpath = tmp % (n, iy)
            newao = yoda.BinnedEstimate1D(ao.xEdges(), newpath)
            for ix in range(1, ao.numBinsX()+1):
                xvis = ao.numBinsX()
                #if ix != xvis:
                #    continue
                b = ao.bin(ix,iy)
                newao.set(ix, b)
            out.append(newao)
        return out
    return ao

