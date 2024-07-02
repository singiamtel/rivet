// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > Omega- Omegabar+
  class BESIII_2022_I2611486 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2611486);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==3334), "UFS");
      // histos
      book(_sigma,1,1,1);
      for (const string& en : _sigma.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // find the Omegas
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p1 : ufs.particles()) {
        bool matched = false;
        // find the children
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        for(const Particle& p2 : ufs.particles()) {
          if(p2.pid() != -p1.pid()) continue;
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
            _sigma->fill(_ecms);
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /femtobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2611486);

}
