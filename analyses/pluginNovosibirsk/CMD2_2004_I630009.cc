// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD2_2004_I630009 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_2004_I630009);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_numPiPiGamma, 1, 1, 1);
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
      // three particles (pi0 pi0 gamma)
      if(ntotal!=3) vetoEvent;
      if(nCount[111]==2 && nCount[22]==1)
	_numPiPiGamma->fill(round(sqrtS()/MeV));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_numPiPiGamma, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _numPiPiGamma;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_2004_I630009);


}
