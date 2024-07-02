// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-pi0
  class TASSO_1986_I228876 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1986_I228876);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,0.9,2.4)) {
        book(_nMeson, "TMP/nMeson", refData(1, 1, 1));
      }
      else {
        throw Error("Invalid CMS energy for TASSO_1986_I228876");
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
      // check the 3 meson final state
      if (ntotal==3 && nCount[PID::PI0]==1 &&
          nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 ) {
        _nMeson->fill(sqrtS()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nMeson, crossSection()/nanobarn/sumOfWeights());
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_nMeson,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nMeson;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1986_I228876);

}
