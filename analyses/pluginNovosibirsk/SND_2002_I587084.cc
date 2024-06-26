// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2002_I587084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2002_I587084);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_numPiPiGamma, "TMP/PiPiGamma", refData(1,1,1));

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
        _numPiPiGamma->fill(sqrtS()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_numPiPiGamma,crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      barchart(_numPiPiGamma,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _numPiPiGamma;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2002_I587084);


}
