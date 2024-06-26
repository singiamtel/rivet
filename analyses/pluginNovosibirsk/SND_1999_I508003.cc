// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- -> pi+pi-pi0
  class SND_1999_I508003 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_1999_I508003);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      book(_num3pi, 1, 1, 1);

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
	_num3pi->fill(round(sqrtS()/MeV));

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_num3pi, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _num3pi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_1999_I508003);


}
