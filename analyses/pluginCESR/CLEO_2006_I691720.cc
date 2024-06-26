// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- -> hadrons
  class CLEO_2006_I691720 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2006_I691720);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=1;ix<18;++ix) {
        if(ix==6||ix==16) continue;
	book(_nMeson[ix], 1,1,ix);
      }
      for (const string& en : _nMeson[1].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      if(ntotal==3) {
	if(nCount[211]==1 && nCount[-211]==1 && nCount[111]==1)
	  _nMeson[1]->fill(_ecms);
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
       	if(p1.children().empty()) continue;
       	if(p1.pid()!=113 && abs(p1.pid())!=313) continue;
     	map<long,int> nRes = nCount;
     	int ncount = ntotal;
     	findChildren(p1,nRes,ncount);
      	if(p1.pid()==113 || p1.pid()==223 || p1.pid()==10113) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==111 && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    if(p1.pid()==113) {
	      _nMeson[2]->fill(_ecms);
	      _nMeson[3]->fill(_ecms);
	    }
	    else if(p1.pid()==223) {
	      _nMeson[5]->fill(_ecms);
	    }
	    else if(p1.pid()==10113) {
	      _nMeson[15]->fill(_ecms);
	    }
	  }
	}
      	else if(p1.pid()==abs(113) || abs(p1.pid())==10113) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  int ipi = p1.pid()>0 ? -211 : 211;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==ipi && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    if(p1.pid()==abs(113)) {
	      _nMeson[2]->fill(_ecms);
	      _nMeson[4]->fill(_ecms);
	    }
	    else {
	      _nMeson[15]->fill(_ecms);
	      _nMeson[17]->fill(_ecms);
	    }
	  }
	}
      	else if(p1.pid()==abs(323)) {
	  if(ncount!=1) continue;
	  bool matched = true;
	  int iK = p1.pid()==323 ? -321 : 321;
	  for(auto const & val : nRes) {
	    if(abs(val.first)==iK && val.second!=1) {
	      matched = false;
	      break;
	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
	  }
       	  if(matched) {
	    if(_ecms=="3.671"s) _nMeson[14]->fill(_ecms);
	  }
	}
	// second unstable particle
	if(abs(p1.pid())!=313 && p1.pid()!=113 && p1.pid()!=223 &&
	   p1.pid()!=333 && p1.pid()!=221 && p1.pid()!=331)
	  continue;
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  const Particle& p2 = ufs.particles()[iy];
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount!=0) continue;
	  bool matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(!matched) continue;
	  if( (p1.pid()==113 && p2.pid()==221) ||
	      (p2.pid()==113 && p1.pid()==221) )
	    _nMeson[7]->fill(_ecms);
	  else if( (p1.pid()==223 && p2.pid()==221) ||
		   (p2.pid()==223 && p1.pid()==221) )
	    _nMeson[8]->fill(_ecms);
	  else if( (p1.pid()==333 && p2.pid()==221) ||
		   (p2.pid()==333 && p1.pid()==221) )
	    _nMeson[9]->fill(_ecms);
	  else  if( (p1.pid()==113 && p2.pid()==331) ||
		    (p2.pid()==113 && p1.pid()==331) )
	    _nMeson[10]->fill(_ecms);
	  else  if( (p1.pid()==223 && p2.pid()==331) ||
		    (p2.pid()==223 && p1.pid()==331) ) {
	    if(_ecms=="3.773"s)_nMeson[11]->fill(_ecms);
          }
	  else  if( (p1.pid()==333 && p2.pid()==331) ||
		    (p2.pid()==333 && p1.pid()==331) ) {
	    if(_ecms=="3.773"s)_nMeson[12]->fill(_ecms);
          }
	  else  if( (p1.pid()==313 && p2.pid()==-313) ||
		    (p2.pid()==313 && p1.pid()==-313) )
	    _nMeson[13]->fill(_ecms);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<18;++ix) {
        if(ix==6||ix==16) continue;
        scale(_nMeson[ix],crossSection()/ sumOfWeights() /picobarn);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nMeson[18];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2006_I691720);


}
