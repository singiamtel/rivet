import yoda, math

def patch(path, ao):
    if "NUSEA_2003_I613362" in path:
        # really only normal hist with mapped axis defined
        newao = yoda.BinnedEstimate1D(ao.xEdges(), ao.path())
        N = newao.numBins()
        # the x-axis has at least one gap -> mask bins!
        mask = [ N - 3 ]
        if any(d in path for d in [ 'd01', 'd02', 'd03', 'd15', 'd16', 'd17', 'd18' ]):
            mask = [ N - 1 ]
        elif any(d in path for d in [ 'd11', 'd24', 'd25', 'd26', 'd29' ]):
            mask = [ N - 4 ]
        elif any(d in path for d in [ 'd19', 'd28' ]):
            mask = [ N - 2 ]
        elif 'd15' in path:
            mask = [ N - 4, N - 1 ]
        newao.maskBin(mask)
        xvals = [ b.xMid() for b in newao.bins() ] # surviving x-values
        yvals = ao.yEdges()
        zvals = ao.zEdges()
        idx = 0
        for x in xvals:
          idx += 1
          if newao.bin(idx).isMasked():
            idx += 1
          # check y-z plane for non-NaN bin value
          skip = False
          for y in yvals:
              for z in zvals:
                  b = ao.binAt(x, y, z)
                  if math.isnan(b.val()):
                      continue
                  newao.set(idx, b)
                  skip = True
                  break # skip rest of z
              if skip:
                  break # skip rest y
        ao = newao
    return ao
