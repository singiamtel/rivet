// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD3_2016_I1395968 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2016_I1395968);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_cKpKmpippim, "TMP/KpKmpippim");
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
      if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1)
	_cKpKmpippim->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sigma = _cKpKmpippim->val();
      double error = _cKpKmpippim->err();
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
    CounterPtr _cKpKmpippim;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMD3_2016_I1395968);


}
