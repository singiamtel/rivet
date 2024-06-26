import yoda
DELPHI_1999_499183_BDIFF=[]
iFilled=0
def patch(path, ao):
    global DELPHI_1999_499183_BDIFF, iFilled
    # fix bin widths
    if "/DELPHI_1999_I499183" and "d39" in path:
        bins=[0.,0.02,0.04,0.06,0.08,0.1,0.14,0.18,0.24,0.34]
        if len(DELPHI_1999_499183_BDIFF)==0 :
            for i in range(1,4) :
                newHist = yoda.BinnedEstimate1D(bins,"/REF/DELPHI_1999_I499183/d39-x01-y0%s"%i)
                newHist.setAnnotation("IsRef",1)
                DELPHI_1999_499183_BDIFF.append(newHist)
        iy = int(ao.path().split("-y")[-1])
        if iy>=1 :
            iFilled +=1
            for ix in range(0,3) :
                DELPHI_1999_499183_BDIFF[ix].set(iy,ao.bin(ix+1))
        if iFilled==9 : return DELPHI_1999_499183_BDIFF
        else : return None
    return ao


