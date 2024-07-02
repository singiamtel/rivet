// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> KKpi
  class CELLO_1989_I266414 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1989_I266414);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if(inRange(sqrtS()/GeV,1.4,4.2)) {
       	book(_nKKPi,"TMP/nKKPi");
      }
      else {
        throw Error("Invalid CMS energy for CELLO_1989_I266414");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if (ntotal!=3) vetoEvent;
      if (nCount[PID::K0S]==1 &&
          ( (nCount[PID::PIPLUS ]==1 && nCount[PID::KMINUS]==1 ) ||
            (nCount[PID::PIMINUS]==1 && nCount[PID::KPLUS]==1 ))) {
        _nKKPi->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nKKPi, crossSection()/nanobarn/sumOfWeights());
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(_nKKPi->val(), _nKKPi->err());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nKKPi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CELLO_1989_I266414);

}
