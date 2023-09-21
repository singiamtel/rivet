// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- -> 3pi+3pi-, 2pi+pi-2pi0, 2pi+pi-K+K-
  class BABAR_2006_I709730 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I709730);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_num3pip3pim,      "TMP/num3pip3pim"     );
      book(_num2pip2pim2pi0,  "TMP/num2pip2pim2pi0" );
      book(_num2pip2pim2KpKm, "TMP/num2pip2pim2KpKm");
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

      if(ntotal!=6) vetoEvent;
      if(nCount[-211]==3 && nCount[211]==3)
	_num3pip3pim->fill();
      else if(nCount[-211]==2 && nCount[211]==2 && nCount[111]==2)
	_num2pip2pim2pi0->fill();
      else if(nCount[-211]==2 && nCount[211]==2 && nCount[321]==1 && nCount[-321]==1)
	_num2pip2pim2KpKm->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ix=1; ix<4; ++ix) {
        double sigma = 0., error = 0.;
        if(ix==1) {
          sigma = _num3pip3pim->val();
          error = _num3pip3pim->err();
        }
        else if(ix==2) {
          sigma = _num2pip2pim2pi0->val();
          error = _num2pip2pim2pi0->err();
        }
        else if(ix==3) {
          sigma = _num2pip2pim2KpKm->val();
          error = _num2pip2pim2KpKm->err();
        }
        sigma *= crossSection()/ sumOfWeights() /nanobarn;
        error *= crossSection()/ sumOfWeights() /nanobarn;
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            mult->bin(1).set(sigma, error);
          }
        }
      }
    }

    /// @}

    // just count the number of events of the types we're looking for
    CounterPtr _num3pip3pim,_num2pip2pim2pi0,_num2pip2pim2KpKm;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I709730);


}
