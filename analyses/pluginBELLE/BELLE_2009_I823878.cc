// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2009_I823878 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I823878);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_nMeson[1], 1, 1, 1);
      book(_nMeson[2], 1, 1, 2);
      book(_nMeson[3], 1, 1, 3);
      book(_nMeson[4], 1, 1, 4);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
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
        if (p.children().empty()) continue;
        if (p.pid()!=333 && p.pid()!=113) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        for (const Particle& p2 : ufs.particles()) {
          if (p2.pid()!=221 && p2.pid()!=331 ) continue;
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
            unsigned int iloc=1;
            if(p.pid()==113 ) iloc+=2;
            if(p2.pid()==331) iloc+=1;
            _nMeson[iloc]->fill(Ecm);
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/ sumOfWeights() /femtobarn;
      scale(_nMeson, sf);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nMeson[5];
    const string Ecm = "10.58";
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BELLE_2009_I823878);


}
