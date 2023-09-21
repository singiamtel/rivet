// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class JADE_1979_I142874 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1979_I142874);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_nHadrons, "TMP/hadrons");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _nHadrons->fill(fs.particles().size());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nHadrons->val()/sumOfWeights();
      double error = _nHadrons->err()/sumOfWeights();
      Estimate1DPtr mult;
      book(mult, 2, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nHadrons;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(JADE_1979_I142874);


}
