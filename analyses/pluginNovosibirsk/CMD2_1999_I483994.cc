// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief e+e- -> 4 pions
  class CMD2_1999_I483994 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_1999_I483994);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_ncharged, 1, 1, 1);
      book(_nneutral, 2, 1, 1);
      book(_nomega  , 3, 1, 1);
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
      if(ntotal==4) {
	if( nCount[211] == 2 && nCount[-211] == 2 )
	  _ncharged->fill(round(sqrtS()/MeV));
	else if( nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 2)
	  _nneutral->fill(round(sqrtS()/MeV));
      }


      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// find the omega
	if(p.pid()==223) {
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // omega pi+pi-
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111) {
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
	    _nomega->fill(round(sqrtS()/MeV));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_ncharged,fact);
      scale(_nneutral,fact);
      scale(_nomega  ,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _ncharged,_nneutral,_nomega;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_1999_I483994);


}
