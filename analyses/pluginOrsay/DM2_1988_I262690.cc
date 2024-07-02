// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- -> K+K-
  class DM2_1988_I262690 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM2_1988_I262690);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nkaon, "TMP/kaon", refData(1, 1, 1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::KPLUS) vetoEvent;
      }
      _nkaon->fill(sqrtS()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nkaon, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_nkaon,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nkaon;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM2_1988_I262690);


}
