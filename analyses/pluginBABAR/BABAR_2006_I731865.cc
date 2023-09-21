// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I731865 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I731865);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_nPhiEta, 1, 1, 2);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
          nRes[child.pid()]-=1;
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
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        if(p.pid()!=333) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        for (const Particle& p2 : ufs.particles()) {
          if (p2.pid()!=221) continue;
          if (p2.parents()[0].isSame(p)) continue;
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if (ncountB!=0) continue;
          bool matched2 = true;
          for (const auto& val : nResB) {
            if (val.second!=0) {
              matched2 = false;
              break;
            }
          }
          if (matched2) {
            _nPhiEta->fill(Ecm);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nPhiEta, crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nPhiEta;
    const string Ecm = "10.58";
    /// @}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I731865);

}
