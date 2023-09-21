def patch(path, ao):
    if "L3_1992_I334954" in path :
        if "d16" in path:
            # remove average bin
            ao.maskBin(ao.indexAt('> 6.0'))
    return ao
