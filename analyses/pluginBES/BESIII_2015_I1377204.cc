// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to J/\psi \pi^0\pi^0$ at energies between 4.19 and 4.42 GeV
  class BESIII_2015_I1377204 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1377204);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nJPsi, "TMP/jpsi");
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
        // find the j/psi
        if(p.pid()==443) {
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          // J/psi pi0pi0
          if(ncount!=2) continue;
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==111) {
              if(val.second !=2) {
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
            _nJPsi->fill();
            break;
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nJPsi->val();
      double error = _nJPsi->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      Estimate1DPtr  mult;
      book(mult, 1, 1, 12);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nJPsi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1377204);


}
