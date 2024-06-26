import yoda
def patch(path, ao):
    if( "TASSO_1988_I263859" in path and "d06" in path) :
        bins=[2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]
        newAO = yoda.BinnedEstimate1D(bins,ao.path())
        for i in range(0,len(bins)) :
            newAO.set(i+1,ao.bin(i+1))
        ao=newAO
    return ao
