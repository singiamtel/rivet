// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> omega pi0pi0
  class BESIII_2021_I1999208 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1999208);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numOmegaPiPi , "TMP/OmegaPiPi");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
	if(child.children().empty()) {
	  --nRes[child.pid()];
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==223)) {
     	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;	    
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==111) {
	    if(val.second !=2) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched) {
	  _numOmegaPiPi->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _numOmegaPiPi->val()* crossSection()/ sumOfWeights() /picobarn;
      double error = _numOmegaPiPi->err()* crossSection()/ sumOfWeights() /picobarn;
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
    CounterPtr _numOmegaPiPi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1999208);

}
