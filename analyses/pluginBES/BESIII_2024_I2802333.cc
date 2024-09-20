// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > K- Xibar+ Lambda0/Sigma0
  class BESIII_2024_I2802333 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2802333);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==3122 or Cuts::abspid==3212 or Cuts::abspid==3312), "UFS");
      // counter
      for(unsigned int ix=0;ix<2;++ix)
        book(_sigma[ix], 1+ix, 1, 1);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int &ncount) {
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
      // first for the Xi
      for (const Particle& xi : ufs.particles(Cuts::abspid==3312)) {
        bool matched = false;
        map<long,int> nRes1=nCount;
        int ncount1 = ntotal;
        findChildren(xi,nRes1,ncount1);
        int sign = xi.pid()/xi.abspid();
        for (const Particle& lam : ufs.particles(Cuts::pid==-sign*3122 or Cuts::pid==-sign*3212)) {
          map<long,int> nRes2=nRes1;
          int ncount2 = ncount1;
          findChildren(lam,nRes2,ncount2);
          if (ncount2!=1) continue;
          matched = true;
          for (auto const& val : nRes2) {
            if (val.first== sign*321) {
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
          if(matched) {
            if(lam.abspid()==3122) _sigma[0]->fill(_ecms);
            else                   _sigma[1]->fill(_ecms);
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor of 1/2 from considering botyh charge states
      scale(_sigma, 0.5*crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2802333);

}
