// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2018_I1637194 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2018_I1637194);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKKpi,  "/TMP/nKKpi" );
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

      if(ntotal==3 && nCount[130]==1 && nCount[310]==1 && nCount[111] ==1)
        _nKKpi->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nKKpi->val();
      double error = _nKKpi->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn;
      Estimate1DPtr mult;
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
    CounterPtr _nKKpi;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2018_I1637194);


}
