// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> phi lambda lambda bar
  class BESIII_2021_I1859248 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1859248);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counter
      book(_c_phiLL, "TMP/c_phiLL");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over phi mesons
      for (const Particle& phi : ufs.particles(Cuts::pid==333)) {
	bool matched = false;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(phi,nRes,ncount);
	// then Lambda baryons
	for (const Particle& lambda : ufs.particles(Cuts::pid==3122)) {
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(lambda,nResB,ncountB);
	  for (const Particle& lambar : ufs.particles(Cuts::pid==-3122)) {
	    map<long,int> nResC = nResB;
	    int ncountC = ncountB;
	    findChildren(lambar,nResC,ncountC);
	    for(auto const & val : nResC) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      _c_phiLL->fill();
	      break;
	    }
	  }
	  if(matched) break;
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sigma = _c_phiLL->val()*fact;
      double error = _c_phiLL->err()*fact;
      Estimate1DPtr  mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_phiLL;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1859248);

}
