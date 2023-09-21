// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > omega pi0 and omega eta
  class BESIII_2020_I1817739 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1817739);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_numOmegaPi, "TMP/OmegaPi");
      book(_nOmegaEta, "/TMP/nOmegaEta");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else
          findChildren(child,nRes,ncount);
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
      for (const Particle& p : ufs.particles(Cuts::pid==223)) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // check for omega pi0
        bool matched = true;
        for(auto const & val : nRes) {
          if(val.first==111) {
            if(val.second!=1) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if(matched) {
          _numOmegaPi->fill();
          break;
        }
        // now for omega eta
        for (const Particle& p2 : ufs.particles(Cuts::pid==221)) {
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if(ncountB!=0) continue;
          matched = true;
          for(auto const & val : nResB) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            _nOmegaEta->fill();
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=1;ix<3;++ix) {
        double sigma,error;
        if(ix==1) {
          sigma = _nOmegaEta->val();
          error = _nOmegaEta->err();
        }
        else {
          sigma = _numOmegaPi->val();
          error = _numOmegaPi->err();
        }
        sigma *= crossSection()/ sumOfWeights() /picobarn;
        error *= crossSection()/ sumOfWeights() /picobarn;
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nOmegaEta;
    CounterPtr _numOmegaPi;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1817739);

}
