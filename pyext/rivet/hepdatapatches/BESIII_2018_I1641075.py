import yoda,math
BESIII_2018_I1641075_h1=None
BESIII_2018_I1641075_h2=None
BESIII_2018_I1641075_h3=None
# subract background from observed correct with efficiency and normalise
def patch(path, ao):
    global BESIII_2018_I1641075_h1, BESIII_2018_I1641075_h2, BESIII_2018_I1641075_h3
    if "BESIII_2018_I1641075" in path:
        if path == "/REF/BESIII_2018_I1641075/d01-x01-y01" : BESIII_2018_I1641075_h1=ao
        if path == "/REF/BESIII_2018_I1641075/d01-x01-y02" : BESIII_2018_I1641075_h2=ao
        if path == "/REF/BESIII_2018_I1641075/d01-x01-y03" : BESIII_2018_I1641075_h3=ao
        if (BESIII_2018_I1641075_h1 and
            BESIII_2018_I1641075_h2 and  BESIII_2018_I1641075_h3) :
            path = "/REF/BESIII_2018_I1641075/d01-x01-y05"
            edges = []
            for val in BESIII_2018_I1641075_h1.xEdges() : edges.append(float(val)-0.0025)
            edges.append(float(BESIII_2018_I1641075_h1.xEdges()[-1])+0.0025)
            newHist = yoda.BinnedEstimate1D(edges,path)
            
            total=0.
            for i in range(0,BESIII_2018_I1641075_h1.numBins()) :
                NS = BESIII_2018_I1641075_h1.bin(i+1).val()
                NB = BESIII_2018_I1641075_h2.bin(i+1).val()
                eff = BESIII_2018_I1641075_h3.bin(i+1).val()
                yVal  = ( NS - NB ) /eff
                err   = math.sqrt( NS+NB)/eff
                newHist.bin(i+1).set(yVal,err)
                total += yVal
            norm = total*0.005
            newHist.scale(1./norm)
            return[ao,newHist]
    return ao
