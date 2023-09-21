// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I700020 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I700020);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nproton, "TMP/proton");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PROTON) vetoEvent;
      }
      _nproton->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sigma = _nproton->val();
      double error = _nproton->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=1;ix<3;++ix) {
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nproton;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I700020);


}
