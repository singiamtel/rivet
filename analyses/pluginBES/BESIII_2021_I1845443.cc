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
      book(_nEta,   1,1,1);
      book(_nOmega, 2,1,1);
      for (const string& en : _nOmega.binning().edges<0>()) {
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
	    _nEta->fill(_ecms);
	  else
	    _nOmega->fill(_ecms);
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      scale(_nEta  ,fact);
      scale(_nOmega,fact);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _nEta,_nOmega;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1845443);

}
