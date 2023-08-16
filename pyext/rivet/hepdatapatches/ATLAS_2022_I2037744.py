
def patch(path, ao):
    if '/REF/ATLAS_2022_I2037744/d36' in path or '/REF/ATLAS_2022_I2037744/d38' in path:
      for i in range(ao.numPoints()):
          if i < 3:
            p.setX(i)
            p.setXErrs((0.5, 0.5))
          elif i < 4:
            p.setX(i+0.5)
            p.setXErrs((1.0, 1.0))
          else:
            p.setX(i+1.5)
            p.setXErrs((1.0, 1.0))
    return ao

