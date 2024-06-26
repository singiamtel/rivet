// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > psi(2S) pi+pi-
  class BESIII_2021_I1933191 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1933191);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nPsi2s, 1, 1, 3);
      for (const string& en : _nPsi2s.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
	// find the psi(2S)
	if(p.pid()==100443) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // omega pi+pi-
	  if(ncount!=2) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211) {
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
	  if(matched) {
	    _nPsi2s->fill(_ecms);
	    break;
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nPsi2s, crossSection()/ sumOfWeights() /picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _nPsi2s;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1933191);

}
