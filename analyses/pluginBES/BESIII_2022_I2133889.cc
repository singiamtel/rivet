// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> omega pi+ pi-
  class BESIII_2022_I2133889 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2133889);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_sigma[0], 1, 1, 1);
      for (unsigned int ix=1; ix<7; ++ix) {
        book(_sigma[ix], 2, 1, ix);
        if (ix < 3) {
          for (const string& en : _sigma[ix-1].binning().edges<0>()) {
            const double end = std::stod(en)*GeV;
            if (isCompatibleWithSqrtS(end)) {
              _ecms[ix-1] = en;
              break;
            }
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int &ncount) {
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
          if (hasOmegaPiPi) _sigma[0]->fill(_ecms[0]);
        }
        // now omega + second resonance
        for (const Particle& p2 : ufs.particles(Cuts::pid==9000221||
                  Cuts::pid==9010221||
                  Cuts::pid==10221||
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
            if      (p2.pid()==9000221 ) _sigma[1]->fill(_ecms[1]);
            else if (p2.pid()==9010221 ) _sigma[2]->fill(_ecms[1]);
            else if (p2.pid()==10221   ) _sigma[3]->fill(_ecms[1]);
            else if (p2.pid()==225     ) _sigma[4]->fill(_ecms[1]);
            break;
          }
        }
        if (hasOmegaPiPi) break;
      }
      if (hasResonance) return;
      // b_1 pi
      for (const Particle& p : ufs.particles(Cuts::abspid==10213)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if(ncount!=1) continue;
        int ipi = -p.pid()/p.abspid()*211;
        hasResonance=true;
        for (const auto& val : nRes) {
          if (val.first==ipi) {
            if (val.second !=1) {
              hasOmegaPiPi = false;
              break;
            }
          }
          else if(val.second!=0) {
            hasResonance = false;
            break;
          }
        }
        if (hasResonance) {
          _sigma[5]->fill(_ecms[1]);
          break;
        }
      }
      if (hasOmegaPiPi && !hasResonance) _sigma[6]->fill(_ecms[1]);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[7];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2133889);

}
