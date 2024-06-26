// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- to 3 pi from BES
  class BESIII_2019_I1773081 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1773081);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      book(_num3pi, "TMP/num3" , refData(1,1,1));

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
      if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==1)
        _num3pi->fill(sqrtS()/GeV);

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_num3pi, crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr  mult;
      book(mult, 1, 1, 1);
      barchart(_num3pi,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _num3pi;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1773081);

}
