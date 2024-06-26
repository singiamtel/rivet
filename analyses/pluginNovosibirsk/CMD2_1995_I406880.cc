// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> hadrons near phi
  class CMD2_1995_I406880 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_1995_I406880);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nKpKm      , "TMP/KpKm"    , refData(1, 1, 1));
      book(_nK0K0      , "TMP/K0K0"    , refData(1, 1, 2));
      book(_n3pi       , "TMP/3pi"     , refData(1, 1, 3));
      book(_numEtaGamma, "TMP/EtaGamma", refData(1, 1, 4));
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
      if(ntotal==2) {
	if(nCount[321]==1 && nCount[-321]==1)
	  _nKpKm->fill(sqrtS()/MeV);
	else if(nCount[130]==1 && nCount[310]==1)
	  _nK0K0->fill(sqrtS()/MeV);
      }
      else if(ntotal==3 && nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 1)
	_n3pi->fill(sqrtS()/MeV);

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
	    _numEtaGamma->fill(sqrtS()/MeV);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      Estimate1DPtr mult;
      scale(_nKpKm,fact);
      book(mult, 1, 1, 1);
      barchart(_nKpKm,mult);
      scale(_nK0K0,fact);
      book(mult, 1, 1, 2);
      barchart(_nK0K0,mult);
      scale(_n3pi,fact);
      book(mult, 1, 1, 3);
      barchart(_n3pi,mult);
      scale(_numEtaGamma,fact);
      book(mult, 1, 1, 4);
      barchart(_numEtaGamma,mult);
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nKpKm,_nK0K0,_n3pi,_numEtaGamma;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_1995_I406880);


}
