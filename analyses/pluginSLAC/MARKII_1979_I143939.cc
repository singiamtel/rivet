// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKII_1979_I143939 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1979_I143939);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      book(_c_DD, "/TMP/sigma_DD");
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // everything else
      else
	_c_hadrons->fill();
      // identified final state with D mesons
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	bool matched = false;
       	const Particle& p1 = ufs.particles()[ix];
       	int id1 = abs(p1.pid());
       	if(id1 != 411 && id1 != 421) continue;
      	// check fs
      	bool fs = true;
      	for (const Particle & child : p1.children()) {
      	  if(child.pid()==p1.pid()) {
      	    fs = false;
      	    break;
      	  }
      	}
      	if(!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
      	// loop over the other fs particles
       	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
       	  const Particle& p2 = ufs.particles()[iy];
       	  fs = true;
       	  for (const Particle & child : p2.children()) {
       	    if(child.pid()==p2.pid()) {
       	      fs = false;
       	      break;
       	    }
       	  }
       	  if(!fs) continue;
       	  if(p2.pid()/abs(p2.pid())==p1.pid()/abs(p1.pid())) continue;
       	  int id2 = abs(p2.pid());
       	  if(id2 != 411 && id2 != 421) continue;
      	  if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
      	    continue;
      	  map<long,int> nRes2 = nRes;
      	  int ncount2 = ncount;
      	  findChildren(p2,nRes2,ncount2);
      	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    _c_DD  ->fill();
	    break;
	  }
	}
	if(matched) break;
      }
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
      book(hadrons, "sigma_hadrons");
      Estimate1DPtr muons;
      book(muons, "sigma_muons"  );
      Estimate1DPtr mult;
      book(mult, 2,1,1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(R.val(), R.errPos());
          hadrons->bin(b.index()).set(sig_h, err_h);
          muons  ->bin(b.index()).set(sig_m, err_m);
        }
      }
      double sigma = _c_DD->val()*fact;
      double error = _c_DD->err()*fact;
      book(mult, 3,1,1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons, _c_DD;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MARKII_1979_I143939);


}
