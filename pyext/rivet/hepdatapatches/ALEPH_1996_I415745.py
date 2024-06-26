import yoda

# removal of average bins and replace bin no with x values
def patch(path,ao) :
    if "ALEPH_1996_I415745" in path:
        if "d01" in path :
            bins=[0.1,0.15,0.2,0.3,0.4,1.]
            newAO = yoda.BinnedEstimate1D(bins,path)
            for i in range(1,len(bins)+1) :
                newAO.set(i,ao.bin(i))
            newAO2 = yoda.BinnedEstimate1D([0.3,1.0],path.replace("x01","x02"))
            newAO2.set(1,ao.bin(ao.numBins()))
            return [newAO,newAO2]
        elif "d02" in path :
            bins=[0.3,0.6,0.9,1.2,1.5]
            newAO = yoda.BinnedEstimate1D(bins,path)
            for i in range(1,len(bins)+1) :
                newAO.set(i,ao.bin(i))
            newAO2 = yoda.BinnedEstimate1D([1.5,10.],path.replace("x01","x02"))
            newAO2.set(1,ao.bin(ao.numBins()-2))
            newAO3 = yoda.BinnedEstimate1D([0.3,10.],path.replace("x01","x03"))
            newAO3.set(1,ao.bin(ao.numBins()-1))
            newAO4 = yoda.BinnedEstimate1D([0.6,10.],path.replace("x01","x04"))
            newAO4.set(1,ao.bin(ao.numBins()))
            return [newAO,newAO2,newAO3,newAO4]
    return ao
