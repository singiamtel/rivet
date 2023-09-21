// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to K^+K^-$ between 2.00 and 3.08 GeV
  class BESIII_2018_I1704558 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1704558);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nkaon, "TMP/kaon");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if(abs(p.pid())!=PID::KPLUS) vetoEvent;
      }
      _nkaon->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nkaon->val();
      double error = _nkaon->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      Estimate1DPtr  mult;
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
    CounterPtr _nkaon;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1704558);


}
