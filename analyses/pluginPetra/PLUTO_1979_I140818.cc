// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PLUTO_1979_I140818 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1979_I140818);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // counters for R
      book(_c_hadrons, "/TMP/sigma_hadrons");
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
      if(nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        // mu+mu- + photons
        vetoEvent;
      }
      else {
        // everything else
        _c_hadrons->fill();
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      Estimate1DPtr hadrons;
      book(hadrons, 2, 1, 1);
      for (auto& b : hadrons->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sig_h, err_h);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1979_I140818);


}
