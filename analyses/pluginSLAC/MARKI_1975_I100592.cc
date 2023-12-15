// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKI_1975_I100592 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1975_I100592);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_nHadrons, "TMP/hadrons");
      book(_nEvent, "TMP/event");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      _nEvent->fill();
      _nHadrons->fill(fs.particles().size());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=1;ix<3;++ix) {
        double sigma,error;
        if(ix==1) {
          sigma = _nEvent->val()*crossSection()/ sumOfWeights() /nanobarn;
          error = _nEvent->err()*crossSection()/ sumOfWeights() /nanobarn;
        }
        else {
          sigma = _nHadrons->val()/sumOfWeights();
          error = _nHadrons->err()/sumOfWeights();
        }
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
    CounterPtr _nHadrons,_nEvent;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKI_1975_I100592);


}
