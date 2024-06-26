// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > omega pi0, rho eta(')
  class BESII_2004_I661567 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2004_I661567);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // counter
      for(unsigned int ix=0;ix<3;++ix)
	book(_sigma[ix],1,1,1+ix);

      for (const string& en : _sigma[0].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int &ncount) {
      for(const Particle& child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // first the omega pi0 final state
      bool matched=false;
      for (const Particle& p : ufs.particles(Cuts::pid==223)) {
      	if(p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p,nRes,ncount);
      	if (ncount==1) {
      	  matched=true;
      	  for (const auto& val : nRes) {
      	    if (val.first==111) {
      	      if (val.second !=1) {
                matched = false;
                break;
      	      }
      	    }
      	    else if (val.second!=0) {
      	      matched = false;
      	      break;
      	    }
      	  }
      	  if(matched) {
      	    _sigma[0]->fill(_ecms);
	    break;
      	  }
      	}
      }
      if (matched) return;
      // rho eta states
      for (const Particle& p : ufs.particles(Cuts::pid==113)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // check for eta/eta'
        for (const Particle& p2 : ufs.particles(Cuts::pid==221 || Cuts::pid==331)) {
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if (ncountB!=0) continue;
          matched = true;
          for (const auto& val : nResB) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
	    if(p2.pid()==221) _sigma[1]->fill(_ecms);
	    else if(_ecms=="3.686") _sigma[2]->fill(_ecms);
            break;
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3; ++ix) scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2004_I661567);

}
