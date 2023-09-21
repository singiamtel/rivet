
def patch(path, ao):
    if "CLEO_1985_I205668" in path and "d01-x01-y02" in path:
        # remove extra point
        ao.maskBin(4)
    return ao
