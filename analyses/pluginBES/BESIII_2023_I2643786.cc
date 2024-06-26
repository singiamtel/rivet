// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> omega pi+ pi-
  class BESIII_2023_I2643786 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2643786);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<7; ++ix) {
        book(_n[ix], "TMP/n_"+toString(ix+1));
      }
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int& ncount) {
      for(const Particle& child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool hasOmegaPiPi=false, hasResonance=false;
      for (const Particle& p : ufs.particles(Cuts::pid==223)) {
      	if (p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
        // first the omega pi pi final state
        if (ncount==2) {
          hasOmegaPiPi=true;
          for (const auto& val : nRes) {
            if (abs(val.first)==211) {
              if (val.second !=1) {
                hasOmegaPiPi = false;
                break;
              }
            }
            else if(val.second!=0) {
              hasOmegaPiPi = false;
              break;
            }
          }
          if (hasOmegaPiPi) _n[0]->fill();
        }
        if (!hasOmegaPiPi) continue;
       	// now omega + second resonance
       	for (const Particle& p2 : ufs.particles(Cuts::pid==9000221 ||
                                                Cuts::pid==9010221 ||
                                                Cuts::pid==10221 ||
                                                Cuts::pid==225)) {
          if (p2.children().empty()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          hasResonance=true;
          if (ncount2!=0) continue;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              hasResonance = false;
              break;
            }
          }
          if (hasResonance) {
            if      (p2.pid()==9000221 ) _n[1]->fill();
            else if (p2.pid()==9010221 ) _n[2]->fill();
            else if (p2.pid()==10221   ) _n[3]->fill();
            else if (p2.pid()==225     ) _n[4]->fill();
            break;
          }
        }
      	if (hasOmegaPiPi) break;
      }
      if (hasResonance || !hasOmegaPiPi) return;
      // b_1 pi
      for (const Particle& p : ufs.particles(Cuts::abspid==10213 || Cuts::abspid==100213)) {
      	if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount!=1) continue;
        int ipi = -p.pid()/p.abspid()*211;
        hasResonance=true;
        for (const auto& val : nRes) {
          if (val.first==ipi) {
            if (val.second !=1) {
              hasOmegaPiPi = false;
              break;
            }
          }
          else if (val.second!=0) {
            hasResonance = false;
            break;
          }
        }
        if (hasResonance) {
          if (p.abspid()==10213) _n[5]->fill();
          else                   _n[6]->fill();
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<7; ++ix) {
        scale(_n[ix], crossSection()/ sumOfWeights() /picobarn);
        unsigned int ih=2, iy=ix;
        if (ix==0) {
          ih=1;
          iy=1;
        }
        Estimate1DPtr mult;
        book(mult, ih, 1, iy);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(_n[ix]->val(), _n[ix]->err());
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _n[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2643786);

}
