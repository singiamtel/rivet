import yoda,math
ARGUS_1997_I420421_hTemp3=None
ARGUS_1997_I420421_hTemp4=None
# add up the partial waves to get the total
def patch(path, ao):
    global ARGUS_1997_I420421_hTemp3,ARGUS_1997_I420421_hTemp4
    if "ARGUS_1997_I420421" in path:
        print(path)
        if "d04" in path :
            if ARGUS_1997_I420421_hTemp4 :
                for i in range(0,ao.numBins()) :
                    newB = ao.bin(i+1)+ARGUS_1997_I420421_hTemp4.bin(i+1)
                    ARGUS_1997_I420421_hTemp4.set(i+1,newB)
                return[ao,ARGUS_1997_I420421_hTemp4]
            else :
                path = ao.path().replace("y01","y03")
                ARGUS_1997_I420421_hTemp4=ao.clone()
                ARGUS_1997_I420421_hTemp4.setPath(path)
        elif "d03" in path :
            if ARGUS_1997_I420421_hTemp3 :
                ioff=0
                if "y04" in path : ioff=13
                for i in range(0,ao.numBins()) :
                    newB = ao.bin(i+1)+ARGUS_1997_I420421_hTemp3.bin(i+ioff+1)
                    ARGUS_1997_I420421_hTemp3.set(i+ioff+1,newB)
                if "y04" in path :
                    return[ao,ARGUS_1997_I420421_hTemp3]
            else :
                path = ao.path().replace("y01","y05")
                ARGUS_1997_I420421_hTemp3=ao.clone()
                ARGUS_1997_I420421_hTemp3.setPath(path)
    return ao
