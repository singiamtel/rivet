// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> K+K-pi0
  class SND_2020_I1806118 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2020_I1806118);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book( _nKKPi, 1, 1, 1);
      book(_nPhiPi,"TMP/nPhiPi", refData(2,1,1));
      for (const string& en : _nKKPi.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // KK pi state
      if(ntotal==3 && nCount[321]==1 &&
	 nCount[-321]==1 && nCount[111]==1)
	_nKKPi->fill(_ecms);
      // phi pi state
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==333)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=1) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(val.first==111) {
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
	  _nPhiPi->fill(sqrtS()/GeV);
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /nanobarn;
      scale(_nKKPi ,fact);
      scale(_nPhiPi,fact);
      Estimate1DPtr  mult;
      book(mult, 2, 1, 1);
      barchart(_nPhiPi,mult);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nKKPi;
    Histo1DPtr _nPhiPi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2020_I1806118);

}
