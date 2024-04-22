import yoda

b_drap_2D    = [ 0., 0.31, 0.62, 0.94, 1.25, 1.56, 1.88, 2.19, 2.50 ]
b_dphi_2D_m  = [ 0., 0.39, 0.79, 1.18, 1.57, 1.96, 2.36, 2.75, 3.14 ]
b_dphi_2D_p1 = [ 0., 1.65, 2.02, 2.40, 2.77, 3.14 ]
b_dphi_2D_p2 = [ 0., 0.31, 0.63, 0.94, 1.26, 1.57, 1.88, 2.20, 2.51, 2.83, 3.14 ]
b_dphi_2D_e  = [ 0., 0.39, 0.79, 1.18, 1.57, 1.96, 2.36, 2.75, 3.14 ]

needs_patching = {
 'd30' : (78, [ b_drap_2D, b_drap_2D,  b_drap_2D,  b_drap_2D, b_drap_2D ]),
 'd33' : (79, [ b_dphi_2D_m, b_dphi_2D_m,  b_dphi_2D_m,  b_dphi_2D_m, b_dphi_2D_m ]),
 'd36' : (80, [ b_dphi_2D_p1, b_dphi_2D_p2, b_dphi_2D_p2 ]),
 'd39' : (81, [ b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e ]),
 'd66' : (82, [ b_drap_2D, b_drap_2D,  b_drap_2D,  b_drap_2D, b_drap_2D ]),
 'd69' : (83, [ b_dphi_2D_m, b_dphi_2D_m,  b_dphi_2D_m,  b_dphi_2D_m, b_dphi_2D_m ]),
 'd72' : (84, [ b_dphi_2D_p1, b_dphi_2D_p2, b_dphi_2D_p2 ]),
 'd75' : (85, [ b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e, b_dphi_2D_e ]),
}

def patch(path, ao):
    hname = path.split('/')[-1]
    ident = hname.split('-')[0]
    if ident in needs_patching and 'y01' in path:
        gidx = 1
        newIdent, bins = needs_patching[ident]
        rtn = [ ao ]
        for subidx, subbins in enumerate(bins):
            newName ='d%i-x01-y0%i' % (newIdent, subidx+1)
            newPath = path.replace(hname, newName)
            newao = yoda.Estimate1D(subbins, newPath)
            for newb in newao.bins():
                oldb = ao.bin(gidx)
                newb.setVal(oldb.val())
                sources = oldb.sources()
                errs = [ (s, oldb.err(s)) for s in sources ]
                for label, epair in errs:
                    newb.setErr(epair, label)
                gidx += 1
            rtn.append(newao)
        return rtn
    return [ ao ]

