// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- > Xi+ Xi-
  class BESIII_2023_I2695411 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2695411);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==3312), "UFS");
      book(_sigmaXi, 1,1,1);

      for (const string& en : _sigmaXi.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p,map<long,int> & nRes, int &ncount) {
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }

      // find the Lambdas and Sigmas
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p1 : ufs.particles()) {
      	bool matched = false;
      	// check fs
      	bool fs = true;
      	for (const Particle& child : p1.children()) {
      	  if (child.pid()==p1.pid()) {
      	    fs = false;
      	    break;
      	  }
      	}
      	if (!fs) continue;
      	// find the children
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
       	for (const Particle & p2 : ufs.particles(Cuts::pid==-p1.pid())) {
      	  // check fs
      	  bool fs = true;
      	  for (const Particle& child : p2.children()) {
      	    if (child.pid()==p2.pid()) {
      	      fs = false;
      	      break;
      	    }
      	  }
      	  if (!fs) continue;
      	  map<long,int> nRes2 = nRes;
      	  int ncount2 = ncount;
      	  findChildren(p2,nRes2,ncount2);
      	  if (ncount2!=0) continue;
      	  matched=true;
      	  for (const auto& val : nRes2) {
      	    if (val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if(matched) {
	    _sigmaXi->fill(_ecms);
      	    break;
      	  }
       	}
       	if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact  = crossSection()/ sumOfWeights() /femtobarn;
      scale(_sigmaXi, fact);
      double mXi = 1.32171;
      double tau  = 0.25*sqr(sqrtS()/mXi);
      double beta = sqrt(1.-1./tau);
      double alpha = 7.2973525693e-3;
      double GeV2fb = 0.3893793721e12;
      double sigma0 = 4.*M_PI*sqr(alpha/sqrtS())*beta*GeV2fb/3.;
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,2);
      for(const auto & b :_sigmaXi->bins()) {
        tmp->bin(b.index()).set(1e3*sqrt(b.sumW()/(sigma0*(1. + 0.5/tau))),
                               1e3*0.5*sqrt(1./b.sumW()/sigma0/(1. + 0.5/tau))*b.errW());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigmaXi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2695411);

}
