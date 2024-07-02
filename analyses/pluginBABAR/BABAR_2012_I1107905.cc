// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi+ pi- J/psi
  class BABAR_2012_I1107905 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1107905);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==443), "UFS");
      book(_sigmaJPsi, "TMP/jpsi", refData(1,1,1));
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
        // J/psi pi+pi-
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
          _sigmaJPsi->fill(sqrtS());
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigmaJPsi,crossSection()/ sumOfWeights() /picobarn);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_sigmaJPsi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaJPsi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1107905);

}
