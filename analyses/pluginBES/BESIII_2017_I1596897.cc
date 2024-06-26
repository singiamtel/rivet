// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > eta h_c
  class BESIII_2017_I1596897 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1596897);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nhc, 1, 1, 1);
      for (const string& en : _nhc.binning().edges<0>()) {
        const double end = std::stod(en)*MeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else
          findChildren(child,nRes,ncount);
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
      // loop over h_c
      for (const Particle& hc : ufs.particles(Cuts::pid==10443)) {
        if(hc.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(hc,nRes,ncount);
        bool matched = false;
        // loop over eta
        for(const Particle & eta : ufs.particles(Cuts::pid==221)) {
          // check eta not child of h_c
          Particle parent=eta;
          while(!parent.parents().empty()) {
            parent=parent.parents()[0];
            if(parent.abspid()==10443) break;
          }
          if(fuzzyEquals(parent.momentum(),hc.momentum())) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(eta,nRes2,ncount2);
          if(ncount2!=0) continue;
          matched=true;
          for(auto const & val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            _nhc->fill(_ecms);
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nhc, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nhc;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1596897);

}
