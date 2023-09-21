import yoda
def patch(path, ao):
    if "MARKII_1985_I207785" in path:
        # really only normal hist with mapped axis defined
        if "d02" in path or "d04" in path :
          newao = yoda.BinnedEstimate1D(ao.xEdges(), ao.path())
          #x = 2.*points[i].x()/29.
          for i in range(1, newao.numBins()+1):
              b = ao.bin(i,i)
              newao.set(i, b)
          ao = newao
    return ao
