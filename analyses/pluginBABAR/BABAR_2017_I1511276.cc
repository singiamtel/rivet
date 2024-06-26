// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> K_S0K_L0 with pi0, eta and 2pi0
  class BABAR_2017_I1511276 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1511276);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=0;ix<4;++ix) {
        book(_sigma[ix], 1+ix, 1, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty() && _ecms[2].empty() && _ecms[3].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      // stable histos
      if( ntotal == 3 && nCount[130] == 1 &&
	  nCount[310] == 1 && nCount[111] == 1) {
	if(!_ecms[0].empty()) _sigma[0]->fill(_ecms[0]);
      }
      else if( ntotal == 4 && nCount[130] == 1 &&
	       nCount[310] == 1 && nCount[111] == 2) {
	if(!_ecms[3].empty()) _sigma[3]->fill(_ecms[3]);
      }
      // unstable particles
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	if(p.pid()!=333 && p.pid()!=221) continue;
	map<long,int> nRes=nCount;
	int ncount = ntotal;
	findChildren(p,nRes,ncount);
	bool matched  = true;
	if(p.pid()==333 && ncount==1) {
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
	    if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
          }
	}
	else if(p.pid()==221 && ncount==2) {
	  for(auto const & val : nRes) {
	    if(val.first==130 || val.first==310) {
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
	    if(!_ecms[2].empty()) _sigma[2]->fill(_ecms[2]);
          }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<4;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[4];
    string _ecms[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1511276);


}
