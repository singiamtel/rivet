// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > Lambda0 Sigmabar0 + c.c.
  class BESIII_2023_I2686032 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2686032);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nLS, 1, 1, 1);
      for (const string& en : _nLS.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int& ncount) {
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
      for (const Particle& p1 : ufs.particles(Cuts::abspid==3122)) {
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
      	findChildren(p1, nRes, ncount);
        int sign = p1.pid()/3122;
       	for (const Particle& p2 : ufs.particles(Cuts::pid==-sign*3212)) {
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
      	  if (matched) {
            _nLS->fill(_ecms);
      	    break;
      	  }
       	}
       	if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nLS, crossSection()/ sumOfWeights() /picobarn);
      const double mL = 1.115683;
      const double mS = 1.192642;
      const double tau  = sqr(sqrtS()/(mL+mS));
      const double beta = sqrt(1.-2.*(sqr(mS)+sqr(mL))/sqr(sqrtS())+sqr((sqr(mS)-sqr(mL))/sqr(sqrtS())));
      const double alpha = 7.2973525693e-3;
      const double GeV2pb = 0.3893793721e9;
      const double sigma0 = 4.*M_PI*sqr(alpha/sqrtS())*beta*GeV2pb;
      BinnedEstimatePtr<string> tmp;
      book(tmp,1,1,2);
      for(const auto & b : _nLS->bins()) {
        tmp->bin(b.index()).set(sqrt(3.*b.sumW()*tau/(sigma0*(0.5 + 0.5*tau))),
                                sqrt(1.5*tau/b.sumW()/sigma0/(1.+2.*tau))*b.errW());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nLS;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2686032);

}
