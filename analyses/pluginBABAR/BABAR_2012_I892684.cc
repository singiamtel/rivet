// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief e+e- > K+K- + pi+pi-, pi0pi0 or K+K-
  class BABAR_2012_I892684 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I892684);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      bool matched = false;
      for(unsigned int ix=0;ix<8;++ix) {
        book(_sigma[ix],1+ix,1,1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            matched=true;
            break;
          }
        }
      }
      if(!matched)
        MSG_ERROR("Beam energy incompatible with analysis.");
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	// K*0
	if(abs(p.pid())==313) {
	  map<long,int> nRes=nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // K* K+/- pi-/+
	  if(ncount !=2 ) continue;
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
	  if(matched==false) continue;
	  if((nCount[321] == 1 && nCount[-321] ==0 &&
	      nCount[211] == 0 && nCount[-211] == 1) ||
	     (nCount[321] == 0 && nCount[-321] ==1 &&
	      nCount[211] == 1 && nCount[-211] == 0)) {
	    if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
          }
	}
	else if(p.pid()==333) {
	  map<long,int> nRes=nCount;
	  int ncount = ntotal;
	  findChildren(p,nRes,ncount);
	  // phi pi+pi-
	  if(ncount==2) {
	    bool matched = true;
	    for(auto const & val : nRes) {
	      if(abs(val.first)==211) {
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
	  for (const Particle& p2 : ufs.particles()) {
	    if(p2.pid()!=9010221&&p2.pid()!=9000221) continue;
	    if(p2.parents()[0].isSame(p)) continue;
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
	      if(p2.pid()==9010221) {
		if(!_ecms[6].empty()) _sigma[6]->fill(_ecms[6]);
		if(!_ecms[3].empty()) _sigma[3]->fill(_ecms[3]);
	      }
	      else {
		if(!_ecms[4].empty()) _sigma[4]->fill(_ecms[4]);
	      }
	    }
	  }
	}
      }
      if(ntotal==4) {
	if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1) {
	  if(!_ecms[0].empty()) _sigma[0]->fill(_ecms[0]);
        }
	else if( nCount[321]==1 && nCount[-321]==1 && nCount[111]==2) {
	  if(!_ecms[5].empty()) _sigma[5]->fill(_ecms[5]);
        }
	else if( nCount[321]==2 && nCount[-321]==2) {
	  if(!_ecms[7].empty()) _sigma[7]->fill(_ecms[7]);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<8;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[8];
    string _ecms[8];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I892684);


}
