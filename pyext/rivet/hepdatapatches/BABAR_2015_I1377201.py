import yoda,math

# removal of average bins and replace bin no with x values
def patch(path, ao):
    if "BABAR_2015_I1377201" in path :
        temp = path.split("/")[-1].split("-")
        ih = int(temp[0][1:])
        iy = int(temp[2][1:])
        output=[ao]
        for ix in range(1,5) :
            newPath = "/REF/BABAR_2015_I1377201/d%02d-x0%s-y0%s" % (ih+6,ix,iy)
            newHist = yoda.BinnedEstimate1D(ao.yEdges(),newPath)
            newHist.setAnnotation("IsRef",1)
            for ip in range(0,4) :
                newHist.set(ip+1,ao.bin(ix,ip+1))
            output.append(newHist)
        return output
    return ao
