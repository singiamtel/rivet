// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2008_I764099 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I764099);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nUps1pipi, 1, 1, 1);
      book(_nUps2pipi, 2, 1, 1);
      book(_nUps3pipi, 3, 1, 1);
      book(_nUps1KK,   4, 1, 1);
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
	if(p.pid() !=   553  &&
	   p.pid() != 100553 &&
	   p.pid() != 200553 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount!=2) continue;
	bool matched = true;
	for(auto const & val : nRes) {
	  if(abs(val.first)==321 || abs(val.first)==211) {
	    continue;
	  }
	  else if(val.second!=0) {
	    matched = false;
	    break;
	  }
	}
	if(matched) {
	  if(nRes[211]==1 && nRes[-211]==1 ) {
	    if(p.pid()==553)
	      _nUps1pipi->fill("10.87"s);
	    if(p.pid()==100553)
	      _nUps2pipi->fill("10.87"s);
	    if(p.pid()==200553)
	      _nUps3pipi->fill("10.87"s);
	  }
	  else if(nRes[321]==1 && nRes[-321]==1) {
	    if(p.pid()==553)
	      _nUps1KK->fill("10.87"s);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      scale(_nUps1pipi,fact);
      scale(_nUps2pipi,fact);
      scale(_nUps3pipi,fact);
      scale(_nUps1KK  ,fact);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nUps1pipi,_nUps2pipi,_nUps3pipi,_nUps1KK;
    /// @}

  };

  RIVET_DECLARE_PLUGIN(BELLE_2008_I764099);

}
