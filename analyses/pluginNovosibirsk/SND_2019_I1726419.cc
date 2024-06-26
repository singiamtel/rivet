// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief SND pi+pi-pi0eta cross section
  class SND_2019_I1726419 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2019_I1726419);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<4;++ix)
        book(_sigma[ix], 1, 1, 1+ix);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms[0] = en;
          if(end>1.6*GeV) _ecms[1]=en;
          break;
        }
      }
      if(_ecms[0].empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool found = false, foundOmegaPhi = false;
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	// eta pi+pi-
	if(ncount==3) {
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==211 || val.first==111 ) {
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
	    _sigma[0]->fill(_ecms[0]);
	    found = true;
	  }
	}
	for (const Particle& p2 : ufs.particles()) {
	  if(p2.pid()!=223 && p2.pid()!=333) continue;
	  map<long,int> nResB = nRes;
	  int ncountB = ncount;
	  findChildren(p2,nResB,ncountB);
	  if(ncountB!=0) continue;
	  bool matched2 = true;
	  for(auto const & val : nResB) {
	    if(val.second!=0) {
	      matched2 = false;
	      break;
	    }
	  }
	  if(matched2) {
	    if(p2.pid()==223)
	      _sigma[1]->fill(_ecms[0]);
	    else if(p2.pid()==333 && !_ecms[1].empty())
	      _sigma[2]->fill(_ecms[1]);
	    foundOmegaPhi=true;
	  }
	}
      }
      if(found && !foundOmegaPhi)
	_sigma[3]->fill(_ecms[0]);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ix=0;ix<4;++ix) {
        scale(_sigma[ix],fact);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[4];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2019_I1726419);


}
