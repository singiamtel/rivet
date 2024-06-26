// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> KS0 K+- pi-+
  class DM1_1982_I176801 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1982_I176801);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKKpi, "/TMP/nKKpi", refData(1,1,1));

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
      if(ntotal!=3) vetoEvent;
      if(nCount[310]==1 &&
	 ((nCount[ 211]==1&&nCount[-321]==1)||
	  (nCount[-211]==1&&nCount[ 321]==1)))
	_nKKpi->fill(sqrtS()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nKKpi, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_nKKpi,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nKKpi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM1_1982_I176801);


}
