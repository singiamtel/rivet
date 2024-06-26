// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief e+ e- > p pbar
  class DM2_1983_I190558 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM2_1983_I190558);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nproton, "TMP/proton", refData(1,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if(fs.particles().size()!=2) vetoEvent;
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=PID::PROTON) vetoEvent;
      }
      _nproton->fill(sqrtS()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nproton, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_nproton,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nproton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM2_1983_I190558);


}
