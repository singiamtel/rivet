// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi- and 3pi0 or 2pi0eta
  class BABAR_2018_I1700745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2018_I1700745);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {


      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      for(unsigned int ix=0;ix<5;++ix)
        book(_sigma[ix], 1+ix, 1, 1);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) {
        if(isCompatibleWithSqrtS(1.075))
          _ecms = "1.075";
        else
          MSG_ERROR("Beam energy incompatible with analysis.");
      }
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
      if(ntotal==5 && nCount[211]==1 && nCount[-211]==1 && nCount[111]==3)
	_sigma[0]->fill(_ecms);
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
     	if(p.children().empty()) continue;
	// find eta/omegas
      	if(p.pid()==221 || p.pid()==223 ) {
       	  map<long,int> nRes = nCount;
       	  int ncount = ntotal;
       	  findChildren(p,nRes,ncount);
	  // eta
	  if(p.pid()==221) {
	    // 2 pi eta
	    if(ncount==2) {
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
		_sigma[1]->fill(_ecms);
	    }
	    // 4 pi eta
	    else if(ncount==4) {
	      bool matched = true;
	      for(auto const & val : nRes) {
		if(abs(val.first)==211) {
		  if(val.second !=1) {
		    matched = false;
		    break;
		  }
		}
		else if(abs(val.first)==111) {
		  if(val.second !=2) {
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
		_sigma[3]->fill(_ecms);
	    }
	    // pi0 omega eta
	    for (const Particle& p2 : ufs.particles()) {
	      if(p2.pid()!=223) continue;
	      map<long,int> nResB = nRes;
	      int ncountB = ncount;
	      findChildren(p2,nResB,ncountB);
	      if(ncountB!=1) continue;
	      bool matched = true;
	      for(auto const & val : nResB) {
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
	      if(matched) {
		_sigma[4]->fill(_ecms);
		break;
	      }
	    }
	  }
	  else {
	    if(ncount!=2) continue;
	    bool matched = true;
	    for(auto const & val : nRes) {
	      if(abs(val.first)==111) {
		if(val.second !=2) {
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
	      _sigma[2]->fill(_ecms);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for (unsigned int ix=0;ix<5;++ix)
        scale(_sigma[ix],fact);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[5];
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BABAR_2018_I1700745);


}
