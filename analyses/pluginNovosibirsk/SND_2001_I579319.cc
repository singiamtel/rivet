// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2001_I579319 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2001_I579319);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npion1, "TMP/pion1");
      book(_npion2, "TMP/pion2");
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
      if(nCount[-211]==2&&nCount[211]==2)
	_npion1->fill();
      else if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==2)
	_npion2->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=1;ix<3;++ix) {
        double sigma,error;
        if(ix==1) {
          sigma = _npion1->val();
          error = _npion1->err();
        }
        else {
          sigma = _npion2->val();
          error = _npion2->err();
        }
        sigma *= crossSection()/ sumOfWeights() /nanobarn;
        error *= crossSection()/ sumOfWeights() /nanobarn;
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _npion1,_npion2;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2001_I579319);


}
