import math

def patch(path, ao):
    if "CELLO_1983_I191415" in path and ("d04" in path or "d05" in path or "d06" in path):
        # calculate average with missing bin from overlap region
        missingVal = 0.016 if 'd04' in path else 0.14 if 'd05' in path else 0.52
        missingErr = 0.01 if 'd04' in path else 0.05 if 'd05' in path else 0.18
        newVal = 0.5*(ao.bin(6).val() + missingVal)
        newErr = math.sqrt(ao.bin(6).errAvg() ** 2 + missingErr ** 2)
        ao.bin(6).set(newVal, newErr)
        if 'd04' not in path:
          missingVal = 0.016 if 'd05' in path else 0.17
          missingErr = 0.016 if 'd05' in path else 0.17
          newVal = 0.5*(ao.bin(7).val() + missingVal)
          newErr = math.sqrt(ao.bin(7).errAvg() ** 2 + missingErr ** 2)
          ao.bin(7).set(newVal, newErr)
    return ao
