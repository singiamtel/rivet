// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> etapi0
  class CRYSTAL_BALL_1986_I217547 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1986_I217547);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,0.683,2.08)) {
        book(_nEtaPi, "TMP/nEta", refData(1, 1, 1));
        if (inRange(sqrtS()/GeV,0.9,1.1)) {
          book(_h_cTheta,2,1,1);
        }
        else if (inRange(sqrtS()/GeV,1.1,1.480)) {
          book(_h_cTheta,2,1,2);
        }
      }
      else {
        throw Error("Invalid CMS energy for CRYSTAL_BALL_1986_I217547");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
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
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==PID::ETA)) {
        if(p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if(ncount !=1 ) continue;
        bool matched = true;
        for(auto const & val : nRes) {
          if(val.first==PID::PI0) {
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
          double cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
          if(cTheta<=0.9)    _nEtaPi->fill(sqrtS()/MeV);
          if(_h_cTheta ) _h_cTheta ->fill(cTheta);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta ,fact);
      scale(_nEtaPi, fact);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_nEtaPi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nEtaPi, _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1986_I217547);

}
