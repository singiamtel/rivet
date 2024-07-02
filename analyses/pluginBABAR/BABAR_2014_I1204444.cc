// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > psi(2S) pi+ pi-
  class BABAR_2014_I1204444 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2014_I1204444);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==100443), "UFS");
      book(_sigmaPsi, "TMP/psi",refData(1,1,1));
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
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // omega pi+pi-
        if (ncount!=2) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
          _sigmaPsi->fill(sqrtS());
          break;
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigmaPsi, crossSection()/ sumOfWeights() /picobarn);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_sigmaPsi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaPsi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2014_I1204444);

}
