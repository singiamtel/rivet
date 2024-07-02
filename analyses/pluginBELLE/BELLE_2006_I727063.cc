// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief gamma gamma -> Lambda Lmabdabar, Sigma0 Sigmabar0
  class BELLE_2006_I727063 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I727063);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==3122 ||	Cuts::abspid==3212), "UFS");
      // histograms
      for(unsigned int ix=0; ix<2; ++ix) {
        book(_sigma[ix],"TMP/c_"+toString(ix+1),refData(1,1,1+ix));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particle Lambda,Lambar;
      // loop over baryon
      bool matched = false;
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty() || p.pid()<0) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
      	// and antibaryons
      	for (const Particle& p2 : ufs.particles()) {
      	  if(p2.children().empty() || p2.pid()>0) continue;
      	  map<long,int> nRes2 = nRes;
      	  int ncount2 = ncount;
      	  findChildren(p2,nRes2,ncount2);
      	  if(ncount2!=0) continue;
      	  matched=2;
      	  for(auto const & val : nRes2) {
      	    if(val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if(matched) {
      	    Lambda=p;
      	    Lambar=p2;
      	    break;
      	  }
      	}
      }
      if(!matched) vetoEvent;
      if(Lambda.pid()!=-Lambar.pid()) vetoEvent;
      double cTheta = abs(Lambda.momentum().z()/Lambda.momentum().p3().mod());
      if(cTheta>0.6) vetoEvent;
      if(Lambda.pid()==3122) _sigma[1]->fill(sqrtS());
      else                   _sigma[0]->fill(sqrtS());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_sigma[ix],fact);
        Estimate1DPtr tmp;
        book(tmp, 1, 1, 1+ix);
        barchart(_sigma[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I727063);

}
