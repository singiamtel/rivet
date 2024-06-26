// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > hadrons at 3.67
  class CLEOC_2005_I654529 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2005_I654529);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      for (unsigned int ix=0; ix<14; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      if (ntotal==3 && nCount[111]==1 && nCount[211]==1 && nCount[-211]==1) {
        _h[0]->fill("3.67"s);
      }
      // unstable particles
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // first the omega, rho, phi  and b_10 pi0 final state
      bool matched=false;
      for (const Particle& p : ufs.particles(Cuts::pid==113 ||
					     Cuts::pid==223 ||
					     Cuts::pid==333 ||
					     Cuts::pid==10113)) {
       	if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount==1) {
          matched=true;
          for (const auto& val : nRes) {
            if (val.first==111) {
              if (val.second !=1) {
                matched = false;
                break;
              }
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            if (p.pid()==113) {
              _h[1]->fill("3.67"s);
              _h[2]->fill("3.67"s);
            }
            else if (p.pid()==223) _h[4]->fill("3.67"s);
            else if (p.pid()==333) _h[5]->fill("3.67"s);
            else if (p.pid()==10113) {
              _h[11]->fill("3.67"s);
              _h[12]->fill("3.67"s);
            }
            break;
          }
        }
        else if (p.pid()!=10113) {
          // check for eta
          for (const Particle& p2 : ufs.particles(Cuts::pid==221)) {
            map<long,int> nResB = nRes;
            int ncountB = ncount;
            findChildren(p2,nResB,ncountB);
            if (ncountB!=0) continue;
            matched = true;
            for (const auto& val : nResB) {
              if (val.second!=0) {
                matched = false;
                break;
              }
            }
            if (matched) {
              if     (p.pid()==113) _h[6]->fill("3.67"s);
              else if(p.pid()==223) _h[7]->fill("3.67"s);
              else if(p.pid()==333) _h[8]->fill("3.67"s);
              break;
            }
          }
        }
      }
      if (matched) return;
      // rho+- b1+- pi-+
      for (const Particle& p : ufs.particles(Cuts::abspid==213 ||
					     Cuts::abspid==10213)) {
       	if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if(ncount!=1) continue;
        int sign = - p.pid()/p.abspid();
        matched=true;
        for (const auto& val : nRes) {
          if (val.first==sign*211 ) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
          if (p.abspid()==213) {
            _h[1]->fill("3.67"s);
            _h[3]->fill("3.67"s);
          }
          else  {
            _h[11]->fill("3.67"s);
            _h[13]->fill("3.67"s);
          }
          break;
        }
      }
      if (matched) return;
      // K* K
      for (const Particle& p : ufs.particles(Cuts::abspid==313 || Cuts::abspid==323)) {
       	if (p.children().empty()) continue;
       	map<long,int> nRes = nCount;
      	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
       	if (ncount!=1) continue;
        int iMeson=0;
        int iStable=0;
       	for (const auto& val : nRes) {
       	  if (val.second==1) {
            iStable+=1;
            iMeson=val.first;
       	  }
          else if (val.second!=0) {
       	    matched = false;
       	    break;
      	  }
      	}
        if (iStable==1) {
          if (p.abspid()==313 || (iMeson==310 || iMeson==130)) {
            _h[9]->fill("3.67"s);
            break;
          }
          else if (p.abspid()==323 and abs(iMeson)==321 && p.pid()*iMeson<0) {
            _h[10]->fill("3.67"s);
            break;
          }
          matched=false;
        }
        else {
          matched=false;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection()/picobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h[14];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2005_I654529);
}
