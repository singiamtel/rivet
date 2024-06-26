// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> K+K-pi0 and KS0 K+-pi-+
  class DM2_1991_I318558 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DM2_1991_I318558);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      for (unsigned int ix=0;ix<2;++ix) {
        book(_nMeson[ix], ix+1, 1, 1);
      }
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
	 _nMeson[0]->fill(round(sqrtS()/MeV));
      else if(nCount[321]==1 &&
	      nCount[-321]==1 && nCount[111]==1)
        _nMeson[1]->fill(round(sqrtS()/MeV));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<2;++ix) {
        scale(_nMeson[ix], crossSection()/ sumOfWeights() /nanobarn);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _nMeson[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DM2_1991_I318558);


}
