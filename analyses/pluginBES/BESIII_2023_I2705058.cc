// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- > phi eta'
  class BESIII_2023_I2705058 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2705058);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");

      // loop over eta mesons
      for (const Particle& p : ufs.particles(Cuts::pid==331)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // loop over phi mesons
        for (const Particle& p2 : ufs.particles(Cuts::pid==333)) {
          if (p2.parents()[0].isSame(p)) continue;
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          bool born    = true;
          bool matched = true;
          for (const auto& val : nResB) {
            if (val.first==22) {
              if (val.second!=ncountB) {
                matched = false;
                break;
              }
              if (val.second!=0) born=false;
            }
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
	    _sigma[0]->fill(_ecms);
	    if(born) _sigma[1]->fill(_ecms);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/picobarn;
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_sigma[ix],fact);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2705058);

}
