// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CRYSTAL_BALL_1986_I238081 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1986_I238081);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_c_D_star, "/TMP/sigma_D_star");
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

      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      bool found = false;
      for (const Particle & p : ufs.particles()) {
	if(abs(p.pid())!=413 && abs(p.pid())!=423) continue;
	bool fs = true;
	for (const Particle & child : p.children()) {
	  if(child.pid()==p.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(fs) {
	  found = true;
	  break;
	}
      }
      if(found)
	_c_D_star->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // R
      Estimate0D R = *_c_hadrons/ *_c_muons;
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Estimate1DPtr hadrons;
      book(hadrons, "sigma_hadrons");
      Estimate1DPtr muons;
      book(muons, "sigma_muons");
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(R.val(), R.errPos());
          hadrons->bin(b.index()).set(sig_h, err_h);
          muons  ->bin(b.index()).set(sig_m, err_m);
        }
      }
      // D*
      fact = crossSection()/ sumOfWeights() /nanobarn;
      double sigma = _c_D_star->val()*fact;
      double error = _c_D_star->err()*fact;
      Estimate1DPtr mult2;
      book(mult2, 2, 1, 1);
      for (auto& b : mult2->bins()) {
      	if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
      	  b.set(sigma, error);
      	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons,_c_D_star;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1986_I238081);


}
