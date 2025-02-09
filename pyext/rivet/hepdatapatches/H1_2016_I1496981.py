
q2axis = yoda.Axis([ 5.5, 8.0, 11.0, 16.0, 22.0, 30.0, 42.0, 60.0, 80.0 ])

def patch(path, ao):
    # fix bin widths
    if 'y01' not in path:  return ao
    ident = int(path.split('/')[-1][1:3])
    if ident > 24 and ident != 51:  return ao
    q = (ident-1) % 8
    for b in ao.bins():
        if ident == 51:
            b.scale(1.0 / b.xWidth() / 2.0)
        else:
            b.scale(1.0 / b.xWidth() / q2axis.width(q+1))
    return ao

