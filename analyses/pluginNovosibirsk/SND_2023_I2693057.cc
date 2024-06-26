// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > omega pi0 -> 4pi
  class SND_2023_I2693057 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2023_I2693057);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==223), "UFS");
      book(_nOmegaPi, 1, 1, 1);
      for (const string& en : _nOmegaPi.binning().edges<0>()) {
        double end = std::stod(en)*MeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // require 4 pi final state
      if (ntotal!=4) vetoEvent;
      if (nCount[211]!=1 || nCount[-211]!=1 || nCount[111]!=2) vetoEvent;
      // omega
      for (const Particle& p : apply<FinalState>(event, "UFS").particles()) {
      	if (p.children().empty()) continue;
       	// find the omega
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
      	// omega pi0
       	if (ncount!=1) continue;
       	bool matched = true;
      	for (const auto& val : nRes) {
      	  if (val.first==111) {
      	    if (val.second !=1) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  else if (val.second!=0) {
      	    matched = false;
      	    break;
      	  }
      	}
      	if (matched) {
      	  _nOmegaPi->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nOmegaPi, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nOmegaPi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2023_I2693057);

}
