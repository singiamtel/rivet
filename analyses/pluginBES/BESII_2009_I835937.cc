// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > light hadrons including KS0
  class BESII_2009_I835937 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2009_I835937);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      if (isCompatibleWithSqrtS(3.773)) {
        book(_sigma[1],1,1,2);
        _ecms="3.773";
      }
      else if (isCompatibleWithSqrtS(3.65)) {
        _ecms = "3.65";
      }
      // histograms
      for (unsigned int ix=0; ix<7; ++ix) {
        if(ix>0 && ix<5) continue;
        book(_sigma[ix], 1, 1, 1+ix);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // first the purely FS particles
      if (ntotal==5 && nCount[310]==1 &&
        ((nCount[-321]==2 && nCount[ 321]==1 && nCount[ 211]==1) ||
				 (nCount[ 321]==2 && nCount[-321]==1 && nCount[-211]==1))) {
        _sigma[0]->fill(_ecms);
      }

      // loop over unstable particles
      for (const Particle& p : ufs.particles(Cuts::pid==221 || Cuts::pid==113 || Cuts::abspid==213)) {
        if (p.children().empty()) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p,nRes,ncount);
        // eta /rho0
        bool matched = false;
        if (p.pid()==221 || p.pid()==113) {
       	  if (ncount==3) {
      	    matched = true;
      	    for (const auto& val : nRes) {
      	      if (abs(val.first)==310) {
                if (val.second!=1) {
                  matched = false;
                  break;
                }
      	      }
              else if(abs(val.first)==211) {
                if (val.second!=1) {
                  matched = false;
                  break;
                }
      	      }
              else if (abs(val.first)==321) {
                if (val.second!=1) {
                  matched = false;
                  break;
                }
      	      }
      	      else if (val.second!=0) {
                matched = false;
                break;
      	      }
            }
            if (matched) {
              if (p.pid()==113) _sigma[6]->fill(_ecms);
              else if (p.pid()==221 &&_sigma[1]) _sigma[1]->fill(_ecms);
            }
          }
        }
        // rho+/-
        else if(p.abspid()==213) {
          if (ncount==2) {
            int sign = p.pid()>0 ? -1 : 1;
            matched=true;
            for (const auto& val : nRes) {
              if (abs(val.first)==310) {
                if (val.second!=1) {
                  matched = false;
                  break;
                }
      	      }
              else if (val.first==321*sign) {
                if (val.second!=1) {
                  matched = false;
                  break;
                }
      	      }
      	      else if(val.second!=0) {
                matched = false;
                break;
      	      }
      	    }
            if (matched) _sigma[5]->fill(_ecms);
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      for (unsigned int ix=0; ix<7; ++ix) {
        if (_sigma[ix]) scale(_sigma[ix],fact);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[7];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2009_I835937);

}
