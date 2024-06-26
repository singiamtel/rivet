// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > eta omega and omega pi0
  class BESIII_2022_I2047667 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2047667);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
        book(_nMeson[ix], 1, 1, 1+ix);
      for (const string& en : _nMeson[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==223)) {
	if(p.children().empty()) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	if(ncount==1) {
	  bool matched=true;
	  for(auto const & val : nRes) {
	    if( val.first==111 ) {
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
	    _nMeson[0]->fill(_ecms);
	    continue;
	  }
	}
	else {
	  // check for eta
	  for (const Particle& p2 : ufs.particles(Cuts::pid==221)) {
	    map<long,int> nResB = nRes;
	    int ncountB = ncount;
	    findChildren(p2,nResB,ncountB);
	    if(ncountB!=0) continue;
	    bool matched = true;
	    for(auto const & val : nResB) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      _nMeson[1]->fill(_ecms);
	      break;
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2; ++ix)
        scale(_nMeson[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nMeson[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2047667);

}
