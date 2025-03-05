import logging
def mergeLatex(dm1,dm2) :
    dm1.latexMode()
    dm2.latexMode()
    newLatex=""
    if len(dm1.latex) != len(dm2.latex) :
        dm1.latex = dm1.latex + r"\mathrm{c.c.}"
        dm2.latex = dm1.latex
        return
    for i in range(0,len(dm1.latex)) :
        if dm1.latex[i]==dm2.latex[i] :
            newLatex+=dm1.latex[i]
        elif dm1.latex[i]=="+" and dm2.latex[i]=="-" :
            newLatex+=r"\pm"
        elif dm1.latex[i]=="-" and dm2.latex[i]=="+" :
            newLatex+=r"\mp"
        elif dm1.latex[i]=="+" and dm2.latex[i]=="0" :
            newLatex+="+/0"
        elif dm1.latex[i]=="0" and dm2.latex[i]=="+" :
            newLatex+="0/+"
        elif dm1.latex[i]=="-" and dm2.latex[i]=="0" :
            newLatex+="-/0"
        elif dm1.latex[i]=="0" and dm2.latex[i]=="-" :
            newLatex+="0/-"
        else :
            logging.debug("problem merging latex %s %s" %(dm1.latex,dm2.latex))
    dm1.latex = newLatex
    dm2.latex = newLatex
