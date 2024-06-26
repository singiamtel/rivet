// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B(*) bar{B}(*) exclusive cross section
  class BELLE_2021_I1859137 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1859137);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
        book(_sigma[ix],1,1,1+ix);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
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
      // extract botton hadrons
      Particles bHadrons=apply<FinalState>(event, "UFS").particles(Cuts::abspid==511 or Cuts::abspid==513 or
								   Cuts::abspid==521 or Cuts::abspid==523);
      for(unsigned int ix=0;ix<bHadrons.size();++ix) {
	long pix = bHadrons[ix].parents()[0].abspid();
	if(pix==511 || pix==413 || pix==521 || pix==523) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(bHadrons[ix],nRes,ncount);
	bool matched=false;
	for(unsigned int iy=ix+1;iy<bHadrons.size();++iy) {
	  long piy = bHadrons[ix].parents()[0].abspid();
	  if(piy==511 || piy==413 || piy==521 || piy==523) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(bHadrons[iy],nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(bHadrons[ix].abspid()==511 ||
	       bHadrons[ix].abspid()==521) {
	      if(bHadrons[iy].pid()==-bHadrons[ix].pid())
		_sigma[0]->fill(_ecms);
	      else
		_sigma[1]->fill(_ecms);
	    }
	    else if(bHadrons[ix].abspid()==513 ||
		    bHadrons[ix].abspid()==523) {
	      if(bHadrons[iy].pid()==-bHadrons[ix].pid())
		_sigma[2]->fill(_ecms);
	      else
		_sigma[1]->fill(_ecms);
	    }
	    break;
	  }
	}
	if(matched) break;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<3;++ix)
        scale(_sigma[ix],fact);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1859137);

}
