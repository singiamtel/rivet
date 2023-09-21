// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class GAMMAGAMMA_1979_I141722 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(GAMMAGAMMA_1979_I141722);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_c_charged, "/TMP/Ncharged");
      book(_c_neutral, "/TMP/Nneutral");
      book(_nHadrons, "/TMP/NHadrons");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0),ncharged(0),nneutral(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(PID::isCharged(p.pid()))
	  ncharged += 1;
	else
	  nneutral += 1;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // everything else
      else {
	if(ntotal==2) vetoEvent;
	_c_hadrons->fill();
	_c_charged->fill(ncharged);
	_c_neutral->fill(nneutral);
	_nHadrons->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate0D R = *_c_hadrons/ *_c_muons;
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Estimate1DPtr hadrons;
      book(hadrons, "sigma_hadrons");
      Estimate1DPtr muons;
      book(muons, "sigma_muons"  );
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(R.val(), R.errPos());
          hadrons->bin(b.index()).set(sig_h, err_h);
          muons  ->bin(b.index()).set(sig_m, err_m);
        }
      }
      scale(_c_charged, 1./_nHadrons->sumW());
      scale(_c_neutral, 1./_nHadrons->sumW());
      for (unsigned int iy=1; iy<3;++iy) {
        double aver(0.),error(0.);
        if(iy==1) {
          aver  = _c_charged->val();
          error = _c_charged->err();
        }
        else {
          aver  = _c_neutral->val();
          error = _c_neutral->err();
        }
        Estimate1DPtr mult;
        book(mult, 2, 1, iy);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(aver, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons,_c_neutral,_c_charged;
    CounterPtr _nHadrons;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(GAMMAGAMMA_1979_I141722);


}
