// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > eta phi
  class BELLE_2022_I2146263 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2146263);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_sigmaPhiEta,"/TMP/_nPhiEta",refData(1,1,1));
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

      // loop over eta mesons
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // loop over phi mesons
        for (const Particle& p2 : ufs.particles(Cuts::pid==333)) {
          if (p2.parents()[0].isSame(p)) continue;
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if (ncountB!=0) continue;
          bool matched = true;
          for (const auto& val : nResB) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) _sigmaPhiEta->fill(sqrtS());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigmaPhiEta,crossSection()/ sumOfWeights()/picobarn);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_sigmaPhiEta,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaPhiEta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2146263);

}
