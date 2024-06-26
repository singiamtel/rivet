// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > KSO KL0 pi0
  class BESIII_2023_I2702520 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2702520);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==313 || Cuts::abspid==315), "UFS");

      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_sigma[ix],1+ix,1,1);
      
      for (const string& en : _sigma[0].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p ,map<long,int> & nRes, int &ncount) {
      for (const Particle& child : p.children()) {
        if (child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
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
      // KSO KL0 pi0 final state
      if( ntotal == 3 && nCount[130] == 1 &&
	  nCount[310] == 1 && nCount[111] == 1)
	_sigma[0]->fill(_ecms);

      // unstable particles
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
       	if (p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
       	bool matched  = true;
        if (ncount!=1) continue;
        unsigned int nk0=0;
        for (const auto& val : nRes) {
          if (val.first==310 || val.first==130) {
            nk0 += val.second;
          }
          else if (val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched && nk0==1) {
	  if(p.abspid()==313) _sigma[1]->fill(_ecms);
	  else                _sigma[2]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	scale(_sigma[ix],crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2702520);

}
