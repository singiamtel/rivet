// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > K0S K0L
  class BABAR_2020_I1769654 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2020_I1769654);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nKSKL    , "TMP/nKSKL");
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

      if(ntotal==2 && nCount[130]==1 && nCount[310]==1)
	_nKSKL->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma =  _nKSKL->val()*crossSection()/ sumOfWeights() /picobarn;
      double error =  _nKSKL->err()*crossSection()/ sumOfWeights() /picobarn;
      Estimate1DPtr  mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nKSKL;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2020_I1769654);

}
