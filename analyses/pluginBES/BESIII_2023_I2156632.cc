// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi+pi- J/psi
  class BESIII_2023_I2156632 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2156632);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)
        book(_sigmaPsi[ix], "TMP/psi-" + toString(ix+1), refData(1,1,1+ix));
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p: fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // find the psis
      for (const Particle& p :  ufs.particles(Cuts::pid==443)) {
        if (p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // psi pi+pi-
        if (ncount!=2) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211) {
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
        if (matched) {
	  for(unsigned int ix=0;ix<3;++ix) _sigmaPsi[ix]->fill(sqrtS()/MeV);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
        scale(_sigmaPsi[ix],crossSection()/ sumOfWeights() /picobarn);
        Estimate1DPtr tmp;
        book(tmp,1,1,1+ix);
        barchart(_sigmaPsi[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaPsi[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2156632);

}
