import yoda

def patch(path, ao):
    if "ALICE_2014_I1244523" in path:
        # 2D hist is really only 1D hist
        # x-axis has gaps and values are in reverse x-order
        if "d17" in path or "d18" in path or "d19" in path or "d20" in path or "d21" in path:
            newao = yoda.BinnedEstimate1D(ao.yEdges(), ao.path())
            xvis = ao.numBinsX()
            for iy in range(1, ao.numBinsY()+1):
                for ix in range(1, ao.numBinsX()+1):
                    if ix != xvis:
                        continue
                    b = ao.bin(ix,iy)
                    newao.set(iy, b)
                    xvis -= 2 # account for gaps
                    break
            ao = newao
    return ao
