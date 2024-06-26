// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > p pbar
  class DM1_1979_I141565 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM1_1979_I141565);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nProton,  "/TMP/nProton", refData(1,1,1) );
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==2 && nCount[2212]==1 && nCount[-2212]==1)
	_nProton->fill(sqrtS()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nProton, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_nProton,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nProton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM1_1979_I141565);


}
