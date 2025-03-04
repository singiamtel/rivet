import re, yoda, logging
try :
    import pdg
except:
    logging.critical("You must have the PDG's python api module installed, see https://pdgapi.lbl.gov/doc/pythonapi.html for details")
    quit()
try :
    from particle import Particle
except:
    logging.critical("Need to install the particle python package, see https://pypi.org/project/particle/")
    quit()
# access to the data
api = pdg.connect()
from .special_cases import *
from .particle_data import *
from .decaymode import *
from .extract_br import *
from .templates import *

# convert pdg strings to pids
def getPID(name) :
    global pdgToNumber
    try :
        num = int(name[0])
        name=name[1:]
    except :
        num = 1
    if name in pdgToNumber :
        return (pdgToNumber[name],num)
    else :
        try :
            part = api.get_particle_by_name(name)
            pid = part.mcid
            pdgToNumber[name] = pid
            return(pid,num)
        except :
            if "(1)" in name :
                return getPID(name.replace("(1)","_1"))
            return(None,None)

# get conjugate pid
def conjugatePID(pid) :
    if pid in latexNames:
        if -pid in latexNames:
            return -pid
        else :
            return pid
    else :
        if Particle.from_pdgid(pid).is_self_conjugate:
            return pid
        else :
            return -pid

def getParticle(pid) :
    return Particle.from_pdgid(pid)

def getPDGParticle(pid) :
    return api.get_particle_by_mcid(pid)


### Class for storage of the decay modes from the PDG
class DecayMode() :

    # initialize from the PDG description
    def __init__(self,iloc=None,desc=None) :
        # variables
        self.imode=iloc
        self.finalState = None
        self.modeType   = None
        self.resonances = []
        self.excluded   = []
        self.subModes   = {}
        self.latex = ""
        if not iloc : return
        # start parsing the decay mode, remove decaying particle name
        i  = desc.find("-->")
        self.incoming = getPID(desc[0:i].strip())[0]
        text = desc[i+3:].strip()
        # various things we can't implement
        if "wave" in text or "I=2" in text or "I=0" in text or "via" in text: return
        # make sure space around arrows
        text=text.replace("-->"," --> ")
        # try and standardize the descrpition
        text=text.replace("no ","non-").replace("non- ","non-")
        text=text.replace("2 ","2").replace("3 ","3").replace("4 ","4").replace("5 ","5")
        text=text.replace("total","")
        text=text.replace("chi_c2","chi_c2 ").replace("chi_c2 (","chi_c2(")
        # determine the type of the mode
        # inclusive modes
        if "anything" in text :
            self.modeType="Inclusive"
            text=text.replace("anything","").strip()
            if " + " in text :
                self.modeType="InclusiveSum"
                text=text.replace(" + ","").strip()
        elif " X" in text and "X(" not in text and "X_" not in text and "Xi" not in text:
            self.modeType="Inclusive"
            text=text.replace(" X","").strip()
        # inclusive semi-leptonic
        elif "semileptonic" in text :
            self.modeType="InclusiveSemiLeptonic"
            text=text.replace("semileptonic","").strip()
        # intermedaite particles in decay
        elif "-->" in text :
            self.modeType="Intermediate"
        # exclusive
        else :
            text=text.replace("(non-resonant)","nonresonant")
            text=text.replace("(non resonant)","nonresonant")
            text=text.replace("(nonresonant)","nonresonant")
            text=text.replace("non-resonant","nonresonant")
            text=text.replace("(NR)"        ,"nonresonant")
            self.modeType="Exclusive"
            if "nonresonant" in text :
                text=text.replace("nonresonant","").strip()
                self.excluded.append(-1)
            elif "(ex." in text :
                i = text.find("(ex.")
                nb=1
                j=i+1
                while nb!=0:
                    if text[j]=="("   : nb+=1
                    elif text[j]==")" : nb-=1
                    j+=1
                inter = text[i+4:j-1].split(",")
                text = text[:i]+text[j+1:]
                for pp in inter:
                    ipid = getPID(pp.strip())
                    if ipid[0] :
                        self.excluded.append(ipid[0])
                        if abs(ipid[0])==311  :
                            self.excluded.append(-ipid[0])
                            self.excluded.append(130)
                            self.excluded.append(310)
        # final clean up
        if text[0]=="(" and text[-1]==")" :
            text=text[1:-1]
        # split up the sub modes
        text = text.split(",")
        particles   = {}
        ## loop over sub decays
        for il,level in enumerate(text) :
            particles2={}
            mode=level.split()
            parentId=0
            if len(mode)>=2 and mode[1]=="-->" :
                (parentId,num) = getPID(mode[0])
                # if fails, try removing the mass from name and check with parent
                if not parentId or parentId not in particles :
                    m2 = mode[0].replace("^","")
                    m3 = m2.replace("()","")
                    m4 = m3.replace("_0*","*0")
                    for ii in range(il-1,-1,-1) :
                        for val in text[ii].split() :
                            temp =re.sub(r"\([0-9]*\)","()",val).replace("^","")
                            temp2=temp.replace("()","")
                            temp3 =re.sub(r"\([0-9]*\)","*",val).replace("^","")
                            if ( m2==temp  or m2==temp [:-1] or
                                 m3==temp2 or m3==temp2[:-1] or
                                 m4==temp2 or m4==temp2[:-1] or
                                 m3==temp3 or m3==temp3[:-1] or
                                 m4==temp3 or m4==temp3[:-1] ) :
                                (parentId,num) = getPID(val)
                                break
                        if(parentId) : break
                if parentId in particles and particles[parentId]==num :
                    self.resonances.append(parentId)
                    del particles[parentId]
                    mode=mode[2:]
                else :
                    return
            resDec=False
            particles3={}
            for name in mode :
                if name == "-->" :
                    resDec=True
                    continue
                (pid,num) = getPID(name)
                if not pid :
                    if "non-" in name :
                        if name=="non-rho" or name=="non-rho(770)":
                            self.excluded.append( 113)
                            self.excluded.append( 213)
                            self.excluded.append(-213)
                            continue
                        else :
                            (pid,num)=getPID(name[4:])
                            if pid :
                                self.excluded.append(pid)
                                continue
                if not pid :
                    logging.debug("get pid fails %s in %s" % (name,desc))
                    return
                if not resDec :
                    if pid in particles2 : particles2[pid]+=num
                    else :                 particles2[pid] =num
                else :
                    if pid in particles3 : particles3[pid]+=num
                    else :                 particles3[pid] =num
            self.subModes[parentId]=particles2
            for key,val in particles2.items() :
                if key in particles :
                    particles[key]+=val
                else :
                    particles[key]=val
            if resDec :
                res = 0
                for key in particles2 :
                    if key in particles3 :
                        particles3[key]-=particles2[key]
                        if particles3[key]==0 : del particles3[key]
                    else :
                        res=key
                if res==0 or particles2[res]!=1 :
                    return
                self.resonances.append(res)
                del particles[res]
                self.subModes[res] = particles3
                for key,val in particles3.items() :
                    if key in particles :
                        particles[key]+=val
                    else :
                        particles[key]=val
        # insert rho decays in needed
        if self.modeType=="Intermediate" :
            if 213 in particles and particles[213]==1:
                del particles[213]
                if 211 in particles :
                    particles[211]+=1
                else :
                    particles[211]=1
                if 111 in particles :
                    particles[111]+=1
                else :
                    particles[111]=1
                self.resonances.append(213)
                self.subModes[213] = {211: 1,111:1}
            if -213 in particles and particles[-213]==1:
                del particles[-213]
                if -211 in particles :
                    particles[-211]+=1
                else :
                    particles[-211]=1
                if 111 in particles :
                    particles[111]+=1
                else :
                    particles[111]=1
                self.resonances.append(-213)
                self.subModes[-213] = {-211: 1,111:1}
            if 113 in particles and particles[113]==1:
                del particles[113]
                if 211 in particles :
                    particles[211]+=1
                else :
                    particles[211]=1
                if -211 in particles :
                    particles[-211]+=1
                else :
                    particles[-211]=1
                self.resonances.append(113)
                self.subModes[113] = {211: 1,-211:1}
        self.finalState=particles

    def FS(self) :
        out=[]
        ncount=0
        for key in sorted(self.finalState) :
            nn = self.finalState[key]
            ncount+=nn
            for i in range(0,nn) :
                out.append(key)
        return (ncount,tuple(out))

    def twoBodyCode(self) :
        if self.FS()[0]!=2 :
            logging.debug("twoBodyCode only generates code for 2 body decays")
            return None
        output=""
        def checkString(iloc,i) :
            if i == 80 :
                return "(meson.children()[%s].pid()==211||meson.children()[%s].pid()==321||meson.children()[%s].pid()==2212)" % (iloc,iloc,iloc)
            elif i == -80:
                return "(meson.children()[%s].pid()==-211||meson.children()[%s].pid()==-321||meson.children()[%s].pid()==-2212)" % (iloc,iloc,iloc)
            else :
                return "meson.children()[%s].pid()==%s" % (iloc,i)


        # standard case, no neutral kaons
        if 311 not in self.finalState and -311 not in self.finalState and 310 not in self.finalState and 130 not in self.finalState:
            ids=[]
            for val in self.finalState :
                ids.append(val)
            if len(self.finalState)==2 :
                output=  "if((%s && %s) ||\n                 " % (checkString(0,ids[0]),checkString(1,ids[1]))+\
                    "(%s && %s)) _h_mult[%s]->fill();\n          "  % (checkString(1,ids[0]),checkString(0,ids[1]),self.imode)
            else :
                output="if(%s && %s) _h_mult[%s]->fill();\n          " %(checkString(0,ids[0]),checkString(1,ids[0]),self.imode)
        else :
            if -311 in self.finalState and 311 in self.finalState :
                output="if(nk0==2) _h_mult[%s]->fill();\n          " %(self.imode)
                return output
            else :
                output="if("
                if   310 in self.finalState : output+="nks0==1&&"
                elif 130 in self.finalState : output+="nkl0==1&&"
                elif -311 in self.finalState or 311 in self.finalState : output+="nk0==1&&"
                for val in self.finalState :
                    if val !=310 and val !=130 and abs(val)!=311 : output+= "(%s || %s) &&" %(checkString(0,val),checkString(1,val))
                output = output[:-2] + ") _h_mult[%s]->fill();\n          " % self.imode
        return output

    def latexSubMode(self,particles,sep) :
        global latexNames
        output=""
        for val,num in particles.items() :
            if num!=1 : output+=str(num)
            try :
                if val not in latexNames :
                    part=Particle.from_pdgid(val)
                    latexNames[val]=part.latex_name
                output+=latexNames[val]+" "
            except :
                logging.debug("Unknown particle ID=%s" % val)
                output+=str(val)+" "
            if val in self.subModes :
                output+="(\\to "+self.latexSubMode(self.subModes[val],sep)+")"
        return output

    def latexMode(self,sep=" ") :
        if self.latex !="" : return self.latex
        output="$%s \\to " % getParticle(self.incoming).latex_name
        if len(self.subModes) <=1 :
            output+=self.latexSubMode(self.finalState,sep)
        else :
            output+=self.latexSubMode(self.subModes[0],sep)
        output+="$"
        self.latex=output
        return output

    def titleString(self) :
        output="Title="
        latex = self.latexMode()
        if self.modeType=="Inclusive" :
            output+="Inclusive %s" % (latex[:-1]+"+X$")
        elif self.modeType=="InclusiveSemiLeptonic" :
            output+="Inclusive %s Semi-leptonic" % latex
        elif self.modeType=="InclusiveSum" :
            latex = self.latexMode("/")
            output+="Inclusive %s" % (latex[:-1]+"+X$")
        elif self.modeType=="Exclusive" or self.modeType=="Intermediate" :
            output+="Exclusive %s " % latex
        else :
            logging.debug("latex fails for mode with %s latex %s" % (mode,latex))
        if len(self.excluded)!=0 :
            if len(self.excluded)==1 and self.excluded[0]==-1 :
                output +="(non-res) "
            elif len(self.excluded)==1 :
                part=Particle.from_pdgid(self.excluded[0])
                output+="(non-$%s$)" %part.latex_name
            elif 113 in self.excluded and 213 in self.excluded and -213 in self.excluded and len(self.excluded)==3 :
                output += "(non-$\\rho$)"
            elif 311 in self.excluded and -311 in self.excluded and 130 in self.excluded and 310 in self.excluded :
                output += " (ex. $K^0"
                for val in self.excluded :
                    if abs(val) ==311 or val==310 or val==130 : continue
                    part=Particle.from_pdgid(val)
                    output+=","+part.latex_name
                output += "$)"
            else :
                logging.debug("latex non fails %s" %self.excluded)
        output+="Branching Ratio"
        return output

    def conjugate(self) :
        ccMode = DecayMode()
        ccMode.imode= self.imode
        ccMode.modeType = self.modeType
        ccMode.finalState={}
        for key,val in self.finalState.items() :
            ccMode.finalState[conjugatePID(key)] = val
        for val in self.resonances :
            ccMode.resonances.append(conjugatePID(val))
        for val in self.excluded :
            ccMode.excluded.append(conjugatePID(val))
        for res,mode in self.subModes.items() :
            if res==0 :
                ccRes=res
            else :
                ccRes = conjugatePID(res)
            newSub={}
            for key,val in mode.items() :
                newSub[conjugatePID(key)] = val
            ccMode.subModes[res] = newSub
        return ccMode

def writeDMCode(n,fs,modes,used,subs,exclusive) :
    initFS=False
    code=""
    UnStable=[]
    for dm in modes :
        # inclusive modes
        if dm.modeType=="Inclusive" :
            if n == 1 :
                if fs[0] in stable :
                    subs["inclusive"]+="_h_mult[%s]->fill(cStable[%s]);\n        " %(dm.imode,fs[0])
                else :
                    subs["inclusive"]+="_h_mult[%s]->fill(cUnStable[%s]);\n        " %(dm.imode,fs[0])
            else :
                subs["inclusive"]+="if ("
                for idcode,num in dm.finalState.items() :
                    if idcode in stable :
                        subs["inclusive"]+="  cStable[%s]>=%s &&" %(idcode,num)
                    else :
                        subs["inclusive"]+="cUnStable[%s]>=%s &&" %(idcode,num)
                subs["inclusive"]=subs["inclusive"][:-3] +(")\n           _h_mult[%s]->fill();\n        "%dm.imode)
        # inclusive sum
        elif dm.modeType == "InclusiveSum" :
            subs["inclusive"]+="_h_mult[%s]->fill("%dm.imode
            if n == 2 and -311 in fs and 311 in fs :
                subs["inclusive"]+="cStable[310]+cStable[130]+"
            else :
                for idcode in fs.items() :
                    if idcode in stable :
                        subs["inclusive"]+="cStable[%s]+" %(idcode)
                    else :
                        subs["inclusive"]+="cUnStable[%s]+" %(idcode)
            subs["inclusive"]=subs["inclusive"][:-1]+");\n        "
        # semi-leptonic
        elif dm.modeType == "InclusiveSemiLeptonic":
            if n != 1 or abs(fs[0])==15:
                logging.debug("problem with inclusive semi-leptonic %s" %fs)
                if dm.imode in used : used.remove(dm.imode)
                continue
            inu = abs(fs[0])+1
            if fs[0]>0 : inu*=-1
            subs["inclusive"] += "if (cStable[%s]>=1 && cStable[%s]>=1)\n           _h_mult[%s]->fill();\n        " %(fs[0],inu,dm.imode)
        else :
            # can't be 1 particle
            if n== 1 :
                logging.debug("problem with mode %s %s %s " %(dm.modeType,dm.imode,fs))
                if dm.imode in used : used.remove(dm.imode)
                continue
            # special for 2 body decays (no intermediates)
            elif n==2 :
                code = dm.twoBodyCode()
                # put kaon modes first as can be multiply matches for same FS
                if code[0:2]=="if" : subs["twobody"] = code+ subs["twobody"]
                else :               subs["twobody"] += code
                continue
            # everything else
            # first particle in loop
            if not initFS:
                for idcode in dm.finalState:
                    if idcode not in stable :
                        UnStable.append(idcode)
                # handle no of neutral kaons first
                code="           if("
                nk0 = 0
                if  311 in dm.finalState : nk0 += dm.finalState[ 311]
                if -311 in dm.finalState : nk0 += dm.finalState[-311]
                if nk0!=0 :
                    code += "(cStable[310]+cStable[130])==%s && " % nk0
                for idcode,num in dm.finalState.items() :
                    if abs(idcode)==311: continue
                    if   idcode== 80 : code += "(cStable[211]+cStable[321]+cStable[2212])==%s && " % num
                    elif idcode==-80 : code += "(cStable[-211]+cStable[-321]+cStable[-2212])==%s && " % num
                    else :             code += "cStable[%s]==%s && " %(idcode,num)
                code = code[:-3]+"){\n"
                initFS = True
            # exclusive mode
            if dm.modeType=="Exclusive" :
                # if excluded particles
                code +="              "
                if len(dm.excluded) !=0 :
                    if len(dm.excluded)==1 and dm.excluded[0]==-1 :
                        code+="if(noIntermediates)"
                    else :
                        code += "if("
                        for ex in dm.excluded :
                            code += "(cUnStable.find(%s)==cUnStable.end() || cUnStable[%s]==0) &&" %(ex,ex)
                        code=code[:-2]+") "
                # fill histo
                code += "_h_mult[%s]->fill();\n" % dm.imode
            # intermediates
            elif  dm.modeType=="Intermediate" :
                code+="              if("
                for res in dm.resonances :
                    code+="cUnStable[%s]==1&&" %res
                if len(dm.excluded) !=0 :
                    for ex in dm.excluded :
                        code += "(cUnStable.find(%s)==cUnStable.end() || cUnStable[%s]==0) &&" %(ex,ex)
                code =code[0:-2]+") _h_mult[%s]->fill();\n" % dm.imode

                pass
    if initFS :
        code +="            }\n"
        if len(UnStable)==0 :
            pid=(0,)
        else :
            pid=tuple(sorted(UnStable))
            if pid not in exclusive :
                exclusive[pid] = {}
        if n not in exclusive[pid] : exclusive[pid][n]=""
        exclusive[pid][n] += code
    return subs


def extractBR(pid,ananame,default_stable) :
    # check we can get the particle we need
    try :
        part = getParticle(pid)
    except :
        logging.critical("Can't find particle with PDG code %s." % pid)
        quit()
    # and extract basic info
    # basic properties for the substitutions
    subs = {}
    name = part.name
    subs["rivetName"] = ananame
    subs["pid"] = pid
    subs["name"] = part.name
    subs["latexName"] = part.latex_name
    # storage of the prongs to make multiplicity distribution
    prongs={}
    # all the decay modes
    all_modes={}
    def addMode(dm) :
        (n,fs) = dm.FS()
        if n not in all_modes : all_modes[n]={}
        if fs not in all_modes[n] : all_modes[n][fs] = [dm]
        else : all_modes[n][fs].append(dm)
    # dict for the yoda hists to output
    aos={}
    # factor to convert to %
    toPercent=100.
    # list of hists which are usable
    used =[]
    used2=[]
    # loop over the decay modes from the PDG
    for bf in getPDGParticle(pid).branching_fractions():
        # convert PDG identifier (doeesn't change with PDG version) to an int
        iloc=int(bf.pdgid.split(".")[1].split("/")[0])
        # skip unaccount used for some particles to get rates to add to 1
        if "Unaccounted" in bf.description : continue
        # get the branching ratio and error (if doesn't exist the skip as limit or no data)
        if bf.value :
            y = bf.value*toPercent
        else :
            continue
        if bf.error :
            dy = bf.error*toPercent
        else :
            if bf.error_negative :
                dy = (-bf.error_negative*toPercent,bf.error_positive*toPercent)
            else :
                continue
        # special for prong (i.e. charged particle no) to convert to distribution
        if "prong" in bf.description :
            i = bf.description.find("prong")-1
            if bf.description[i]==" " or bf.description[i]=="-": i-=1
            n = int(bf.description[i])
            prongs[n] = (y,dy)
            used2.append(iloc)
            continue
        # otherwise create yoda object for BR
        if iloc<100 :
            path = "/REF/%s/d01-x01-y%02d" % (ananame,iloc)
        else :
            path = "/REF/%s/d01-x01-y%03d" % (ananame,iloc)
        newHist = yoda.Estimate0D(path)
        newHist.setAnnotation("IsRef",1)
        newHist.set(y,dy)
        aos[path] = newHist
        # parse the PDG description to give the decay mode
        text = bf.description
        text = text.replace(" **","")
        modes=[]
        # check for charge conjugate
        CC = False
        if "c.c" in text :
            text=text.replace("+c.c.","")
            text=text.replace(" + c.c.","")
            text=text.replace("(or c.c.)","")
            text=text.replace("or c.c.","")
            text=text.replace("+ c.c","")
            if("c.c." in text) :
                logging.debug("Charge Conjuate not handled for mode %s" %bf)
                continue
            CC = True
        # x B() often found in B+ modes
        if "x B(" in text :
            istart = text.find("x B(")+4
            icount=1
            for i in range(istart,len(text)) :
                if   text[i]=="(" : icount+=1
                elif text[i]==")" : icount-=1
                if icount==0 : break
            text2 = text[:istart-4] + " , " + text[istart:i]
            if i+1<len(text) : text2+= text[i+1:]
            text=text2
        # sometimes n() around many particles
        for i in range(2,10) :
            text=text.replace(" %s ("%i," %s(" % i)
            il = text.find(" %s("%i)
            while il>=0 :
                newText=text[:il]
                text=text[il+3:]
                ib=1
                for il in range(0,len(text)) :
                    if   text[il]=="(" : ib+=1
                    elif text[il]==")" : ib-=1
                    if(ib==0) : break
                endText=text[il+1:]
                for val in text[0:il].split() :
                    newText+= " %s%s " % (i,val)
                newText +=endText
                text=newText
                il = text.find(" %s("%i)
        # sometimes multiple charge states
        PlusMinus=False
        if "+-" in text or "+0" in text or "-0" in text :
            if "+-" in text or "-+" in text:
                modes.append(text.replace("+-","+").replace("-+","-"))
                modes.append(text.replace("+-","-").replace("-+","+"))
            else :
                modes.append(text.replace("+0","+").replace("0+","0"))
                modes.append(text.replace("+0","0").replace("0+","+"))
            PlusMinus=True
        else :
            modes.append(text)
        # special cases
        if pid in special_modes:
            if iloc in special_modes[pid]:
                CC = False
                modes=[]
                for val in special_modes[pid][iloc]:
                    modes.append(val)
        # loop over strings
        dms=[]
        worked=True
        # convert to DecayMode objects
        for val in modes :
            dms.append(DecayMode(iloc,val))
            if not dms[-1].finalState : worked=False
        if not worked :
            del aos[path]
            continue
        # merge latex strings if multiply charge states
        if PlusMinus :
            mergeLatex(dms[0],dms[1])
        for dm in dms :
            addMode(dm)
            if pid in special_latex and iloc in special_latex[pid] :
                dm.latex = special_latex[pid][iloc]
        # or CC
        if CC :
            for dm in dms :
                dm_cc = dm.conjugate()
                ll = dm.latexMode()[:-1] + " (+\\text{c.c.})$"
                dm.latex    = ll
                dm_cc.latex = ll
                addMode(dm_cc)
        # append to list of used brs
        used.append(iloc)
    # special for multiplicity dists (D0/tau others?)
    if len(prongs)>0 :
        path = "/REF/%s/d02-x01-y01" % ananame
        newHist = yoda.BinnedEstimate1D(list(prongs.keys()),path)
        newHist.setAnnotation("IsRef",1)
        aos[path] = newHist
        for key,val in prongs.items():
            newHist.binAt(key).set(val[0],val[1])
    # subsitutions for the template files
    # set no of histos
    subs["nHist"] = max(used)+1
    # and ones sucessfully parsed
    subs["histIDs"] = ""
    for i,hid in enumerate(sorted(used)) :
        if i>0 and (i%10)==0 :  subs["histIDs"] +="\n                "
        subs["histIDs"] += str(hid) +","
    subs["histIDs"] = subs["histIDs"][:-1]
    # storage for code for the different types of mode
    subs["inclusive"] =""
    subs["twobody"]=""
    exclusive={ (0,) : {} }
    # write the code
    for n,val in all_modes.items() :
        for fs,modes in val.items() :
            writeDMCode(n,fs,modes,used,subs,exclusive)
    # code for dists if needed
    subs["finalize"]=""
    subs["initialize"]=""
    subs["histos"]=""
    subs["charged"]=""
    if len(prongs) !=0 :
        subs["initialize"]+="book(_h_charged,2,1,1);\n        "
        subs["finalize"]+="// 2 due binning of charged dist\n      normalize(_h_charged,100.);\n        "
        subs["histos"] +="BinnedHistoPtr<int> _h_charged;\n        ";
        subs["charged"] +="// first the charged multiplicity\n        _h_charged->fill(nCharged(meson));\n"
    subs["default_stable"]=",".join(str(x) for x in  default_stable)
    subs["exclusive"]=""
    for i,pid2 in enumerate(sorted(exclusive)) :
        start=True
        # set up in unstable particle
        if pid2 !=(0,) :
            temp=",".join(str(x) for x in default_stable)
            hasK0=False
            for val in pid2 :
                if abs(val)==311 :
                    hasK0=True
                    continue
                if abs(val)==80 : continue
                temp+=",%s" % val
            if hasK0 and 310 not in pid2 and 310 not in default_stable:
                temp+=",310"
            subs["exclusive"]+="        // Modes with %s in final state " % temp
            subs["exclusive"]+="\n        unsigned int nstable%s(0);\n        map<int,unsigned int> cStable%s, cUnStable%s;\n        "%(i,i,i)
            subs["exclusive"]+="findDecayProducts(meson,cStable%s,cUnStable%s,nstable%s,set<int>({%s}));\n        noIntermediates = cUnStable.empty();\n        " % (i,i,i,temp)
        loc=str(i)
        if pid2 ==(0,) : loc=""
        temp=""
        for nstable,code in exclusive[pid2].items() :
            if start :
                temp+="if     (nstable%s==%s) {\n" % (loc,nstable)
                start=False
            else :
                temp+="        else if(nstable%s==%s) {\n" % (loc,nstable)
            temp+=code.replace("cStable","cStable%s"%loc).replace("cUnStable","cUnStable%s"%loc)
            temp+="        }\n"
        if pid2 ==(0,) :
            subs["exclusive"] = temp+subs["exclusive"]
        else :
            subs["exclusive"] += temp
    # titles for the plot file
    titles={}
    for n,val in all_modes.items() :
        for fs,modes in val.items() :
            for dm in modes :
                titles[dm.imode] = dm.titleString()
    # plot file
    plot = ""
    if len(prongs) > 0 :
        plot += "BEGIN PLOT /%s/d02-x01-y01\nTitle=Charged particle multiplicity in $ %s $ decays\n" % (ananame,subs["latexName"])
        plot += "XLabel=$n_{\\mathrm{charged}}$\nYLabel=$P(n_{\\mathrm{charged}})$ [$\\% $]\nLogY=0\nEND PLOT\n"
    plot += ("BEGIN PLOT /%s" %ananame)+"/d01-x01\nXLabel=\nYLabel=BR [\\%]\nLogY=0\nEND PLOT\n"
    for i in sorted(titles) :
        if i<100 :
            plot +="BEGIN PLOT /%s/d01-x01-y%02d\n" %(ananame,i)
        else :
            plot +="BEGIN PLOT /%s/d01-x01-y%03d\n" %(ananame,i)
        plot +=titles[i]
        plot +="\nEND PLOT\n"
    # modes not handled
    logging.debug("The following modes were not implemented")
    for bf in getPDGParticle(pid).branching_fractions():
        iloc=int(bf.pdgid.split(".")[1].split("/")[0])
        if "wave" in bf.description : continue
        if "Unaccounted" in bf.description : continue
        if not bf.value or not bf.error_negative: continue
        if iloc not in used and iloc not in used2:
            logging.debug("Fails PDG Identifier=%s, line in table=%s, description=%s" %(iloc,bf.mode_number, bf))
        else :
            logging.debug("Works PDG Identifier=%s, line in table=%s, description=%s" %(iloc,bf.mode_number, bf))
    # return the results
    return (aos,subs,plot)
