// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- -> pi+pi-
  class SND_2005_I686349 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2005_I686349);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, 1,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill(round(sqrtS()/MeV));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npion,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _npion;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2005_I686349);


}
