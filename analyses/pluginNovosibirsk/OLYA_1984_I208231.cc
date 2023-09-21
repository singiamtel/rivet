// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class OLYA_1984_I208231 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OLYA_1984_I208231);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, "TMP/pion");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _npion->val();
      double error = _npion->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn;
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/MeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _npion;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OLYA_1984_I208231);


}
