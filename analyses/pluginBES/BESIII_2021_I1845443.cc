// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- -> p pbar eta/omega
  class BESIII_2021_I1845443 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1845443);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nEta, "/TMP/nEta" );
      book(_nOmega, "/TMP/nOmega" );
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // find the eta/omega
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle& p : ufs.particles(Cuts::pid==PID::ETA or Cuts::pid==PID::OMEGA)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	bool matched=true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==PID::PROTON) {
	    if(val.second!=1) {
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
	  if(p.pid()==PID::ETA)
	    _nEta->fill();
	  else
	    _nOmega->fill();
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
        double sigma,error;
        if(ix==1) {
          sigma = _nEta->val();
          error = _nEta->err();
        }
        else {
          sigma = _nOmega->val();
          error = _nOmega->err();
        }
        sigma *= crossSection()/ sumOfWeights() /picobarn;
        error *= crossSection()/ sumOfWeights() /picobarn;
        Estimate1DPtr  mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _nEta,_nOmega;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1845443);

}
