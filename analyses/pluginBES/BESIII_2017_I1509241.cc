// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to p \bar{p} \pi^0$ at energies between 4.008 and 4.600 GeV
  class BESIII_2017_I1509241 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1509241);


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
      if(fs.particles().size()!=3) vetoEvent;
      for (const Particle& p :  fs.particles()) {
        if(abs(p.pid())!=PID::PROTON&& p.pid()==PID::PI0) vetoEvent;
      }
      _nproton->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nproton->val();
      double error = _nproton->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      Estimate1D temphisto(refData(1, 1, 1));
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
    CounterPtr _nproton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1509241);

}
