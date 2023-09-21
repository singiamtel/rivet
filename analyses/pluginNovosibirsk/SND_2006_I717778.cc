// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2006_I717778 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2006_I717778);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numEtaGamma, "TMP/EtaGamma");
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
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the omega
	if(p.pid()==221) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
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
	    _numEtaGamma->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double sigma = _numEtaGamma->val();
      double error = _numEtaGamma->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn;
      for (unsigned int ix=1;ix<3;++ix) {
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/MeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}

    /// @name Histograms
    /// @{
    CounterPtr _numEtaGamma;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2006_I717778);


}
