// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for D0 D*- pi+ +c.c
  class BESIII_2018_I1685535 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1685535);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nD0, "/TMP/nD0");
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


      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
        const Particle& p1 = ufs.particles()[ix];
        if(abs(p1.pid())!=421) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        bool matched=false;
        int id2 = p1.pid()>0 ? -413 :  413;
        int ipi = p1.pid()>0 ?  211 : -211;
        for(unsigned int iy=0;iy<ufs.particles().size();++iy) {
          if(ix==iy) continue;
          const Particle& p2 = ufs.particles()[iy];
          if(p2.pid()!=id2) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if(ncount2!=1) continue;
          matched=true;
          for(auto const & val : nRes2) {
            if(val.first==ipi) {
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
          if(matched) break;
        }
        if(matched)
          _nD0->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nD0->val();
      double error = _nD0->err();
      sigma *= crossSection()/ sumOfWeights() /nanobarn;
      error *= crossSection()/ sumOfWeights() /nanobarn;
      for (unsigned int ix=1;ix<3;++ix) {
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
    CounterPtr _nD0;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1685535);

}
