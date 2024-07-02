// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for e+e- > gamma eta
  class SND_2014_I1275333 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2014_I1275333);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numEtaGamma, "TMP/EtaGamma", refData(1,1,1));
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
	if(p.children().empty()) continue;
	// find the eta
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta gamma
	if(ncount!=1) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(val.first==22) {
	    if(val.second !=1) {
	      matched = false;
	      break;
	    }
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched)
	  _numEtaGamma->fill(sqrtS()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_numEtaGamma, crossSection()/ sumOfWeights() /picobarn);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_numEtaGamma,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _numEtaGamma;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2014_I1275333);


}
