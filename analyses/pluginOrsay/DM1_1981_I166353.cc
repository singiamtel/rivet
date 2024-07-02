// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > 3pi+3pi-
  class DM1_1981_I166353 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1981_I166353);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_num3pip3pim, "TMP/num3pip3pim",refData(1,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	if(abs(p.pid())!=211) vetoEvent;
	++ntotal;
      }
      if(ntotal!=6) vetoEvent;
      _num3pip3pim->fill(sqrtS()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_num3pip3pim,crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_num3pip3pim,mult);
    }
    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _num3pip3pim;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(DM1_1981_I166353);


}
