// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD3_2016_I1444990 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2016_I1444990);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      book(_nK0K0, "TMP/K0K0");

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
      if(ntotal==2 &&
	 nCount[130]==1 && nCount[310]==1)
	_nK0K0->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nK0K0->val();
      double error = _nK0K0->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn;
      Estimate1DPtr mult;
      book(mult, 1, 1, 6);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/MeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nK0K0;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2016_I1444990);


}
