// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > K+K-pi+pi-
  class DM1_1982_I169382 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1982_I169382);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_cKpKmpippim, "TMP/KpKmpippim",refData(1,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }

      if(ntotal!=4) vetoEvent;
      if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1)
	_cKpKmpippim->fill(sqrtS()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_cKpKmpippim, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_cKpKmpippim,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _cKpKmpippim;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM1_1982_I169382);


}
