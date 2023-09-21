// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD3_2016_I1385598 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2016_I1385598);


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
      if (fs.particles().size() != 2) vetoEvent;
      for (const Particle& p : fs.particles()) {
        if (p.abspid() != PID::PROTON) vetoEvent;
      }
      _nproton->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nproton, crossSection()/ sumOfWeights() /nanobarn);
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 6);
      for (auto& b : mult->bins()) {
        const double sqrts = std::stod(b.xEdge().substr(0, 4));
        if (isCompatibleWithSqrtS(sqrts*MeV)) {
          b.set(_nproton->val(), _nproton->err());
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
  RIVET_DECLARE_PLUGIN(CMD3_2016_I1385598);


}
