// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief e+ e- > pi+pi-
  class KLOE2_2017_I1634981 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KLOE2_2017_I1634981);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion, "TMP/pion", refData(1, 1, 1));
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PIPLUS) vetoEvent;
      }
      _npion->fill(sqr(sqrtS()/GeV));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _npion, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr  mult;
      book(mult, 1, 1, 1);
      barchart(_npion,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _npion;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(KLOE2_2017_I1634981);


}
