// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > pi+ pi- h_c
  class BESIII_2017_I1494065 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1494065);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nhc, "TMP/h_c");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
          --nRes[child.pid()];
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
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        // find the h_c
        if(p.pid()==10443) {
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          // omega pi+pi-
          if(ncount!=2) continue;
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==211) {
              if(val.second !=1) {
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
            _nhc->fill();
            break;
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nhc->val();
      double error = _nhc->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int  ix=1;ix<3;++ix) {
        Estimate1DPtr  mult;
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
    CounterPtr _nhc;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1494065);

}
