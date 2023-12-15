// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class GAMMAGAMMA_1973_I84794 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(GAMMAGAMMA_1973_I84794);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // everything else
      else
	_c_hadrons->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate0D R = *_c_hadrons/ *_c_muons;
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Estimate1DPtr hadrons;
      book(hadrons, 3,1,1);
      Estimate1DPtr muons;
      book(muons, "sigma_muons"  );
      Estimate1DPtr mult;
      book(mult, 4,1,1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(R.val(), R.errPos());
          hadrons->bin(b.index()).set(sig_h, err_h);
          muons  ->bin(b.index()).set(sig_m, err_m);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(GAMMAGAMMA_1973_I84794);


}
