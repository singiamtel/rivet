// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for D*0 D*- pi+ + c.c.
  class BESIII_2023_I2645388 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2645388);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_sigma[ix], 1+ix, 1, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
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
      for (const Particle & p1 : ufs.particles(Cuts::abspid==423)) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        bool matched=false;
      	int id2 = p1.pid()>0 ? -413 :  413;
      	int ipi = p1.pid()>0 ?  211 : -211;
       	for (const Particle & p2 : ufs.particles(Cuts::pid==id2)) {
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=1) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (val.first==ipi) {
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
      	  if (matched) break;
      	}
      	if (matched) {
          for (unsigned int ix=0; ix<2; ++ix) {
            if (!_ecms[ix].empty()) _sigma[ix]->fill(_ecms[ix]);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2645388);

}
