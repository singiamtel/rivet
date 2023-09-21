// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DM1_1981_I166353 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1981_I166353);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_num3pip3pim, "TMP/num3pip3pim");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=211) vetoEvent;
	++ntotal;
      }
      if(ntotal!=6) vetoEvent;
      _num3pip3pim->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sigma = _num3pip3pim->val()*fact;
      double error = _num3pip3pim->err()*fact;
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }
    /// @}

    /// @name Histograms
    /// @{
    CounterPtr _num3pip3pim;
    /// @}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DM1_1981_I166353);


}
