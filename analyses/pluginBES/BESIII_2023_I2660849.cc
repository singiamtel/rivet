// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > Ds*+ Ds*-
  class BESIII_2023_I2660849 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2660849);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==433), "UFS");
      // Book histograms
      book(_sigma,1,1,1);
      double delta=1e30;
      for (const string& en : _sigma.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          if(abs(end-sqrtS())<delta) {
            _ecms = en;
            delta = abs(end-sqrtS());
          }
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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
      // final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // unstable charm analysis
      Particles DS = apply<UnstableParticles>(event, "UFS").particles();
      for (unsigned int ix=0; ix<DS.size(); ++ix) {
        const Particle& p1 = DS[ix];
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
        bool matched=false;
        // loop over the anti particles
        for (unsigned int iy=ix+1; iy<DS.size(); ++iy) {
          const Particle& p2 = DS[iy];
          if (p2.pid() != -p1.pid()) continue;
          fs = true;
          for (const Particle& child : p2.children()) {
            if (child.pid()==p2.pid()) {
              fs = false;
              break;
            }
          }
          if (!fs) continue;
          if (!p2.parents().empty() && p2.parents()[0].pid()==p1.pid()) {
            continue;
          }
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
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/ sumOfWeights()/picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2660849);

}
