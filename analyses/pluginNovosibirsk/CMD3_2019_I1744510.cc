// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > eta pi+pi-
  class CMD3_2019_I1744510 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2019_I1744510);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numEtaPiPi, 2, 1, 1);
      for (const string& en : _numEtaPiPi.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the omega
	if(p.pid()==221) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // eta pi+pi-
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
	  if(matched)
	    _numEtaPiPi->fill(_ecms);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_numEtaPiPi, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _numEtaPiPi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2019_I1744510);


}
