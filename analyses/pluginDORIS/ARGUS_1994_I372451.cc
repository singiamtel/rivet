// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> phi rho or omega
  class ARGUS_1994_I372451 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1994_I372451);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV, 1.5, 3.5)) {
        book(_nMeson, 1,1,1);
        book(_avg, 2,1,1);
      }
      else {
        throw Error("Invalid CMS energy for ARGUS_1994_I372451");
      }
    }

    void findChildren(const Particle& p, map<long,int>& nRes, int& ncount) {
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles phis = ufs.particles(Cuts::pid==333);
      Particles rhos = ufs.particles(Cuts::pid==113 or Cuts::pid==223);
      for (const Particle& phi : phis) {
        bool matched=false;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(phi,nRes,ncount);
        for (const Particle & rho : rhos) {
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(rho,nRes2,ncount2);
          if (ncount2 !=0 ) continue;
          matched = true;
          for (auto const & val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            if (rho.pid()==113) {
              _nMeson->fill(sqrtS()/GeV);
            }
            else {
              _avg->fill("1.9 - 2.3"s);
            }
            break;
          }
        }
        if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/nanobarn/sumOfWeights();
      scale(_nMeson, sf);
      scale(_avg, sf);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nMeson;
    BinnedHistoPtr<string> _avg;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1994_I372451);

}
