import yoda,math

# removal of average bins and replace bin no with x values
def patch(path, ao):
    if "BABAR_2016_I1391152" in path :
        mB = 5.27966
        mK = 0.89167
        avg = [ 1., 6. ]
        bins = [ 0.1, 2., 4.3, 8.12, 10.11, 12.89, 14.21, (mB-mK)**2 ]
        newao = yoda.BinnedEstimate1D(bins, ao.path())
        newao.maskBin(4) # gap between 8.12 and 10.11
        newao.maskBin(6) # gap between 12.89 and 14.21
        avg = yoda.BinnedEstimate1D(avg, ao.path().replace("x01","x02"))
        for i, b in enumerate(ao.bins()):
            if i == 5:
                newao.set(7, b)
            elif i == 4:
                newao.set(5, b)
            elif i:
                newao.set(i, b)
            else:
                avg.set(i+1, b)
        return [ newao, avg ]
    return ao
