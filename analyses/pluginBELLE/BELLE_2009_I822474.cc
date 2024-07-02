// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> eta pi0
  class BELLE_2009_I822474 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I822474);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,0.84,4.)) {
        book(_sigmaEtaPi,"TMP/nEtaPi",refData(1,1,1));
        double sMin=0.84, step=0.02;
        unsigned int ihist=2;
        while (sMin<4.) {
          if (inRange(sqrtS()/GeV, sMin, sMin+step)) {
            break;
          }
          sMin+=step;
          ihist+=1;
          if (fuzzyEquals(1.6, sMin)) step=0.04;
          else if (fuzzyEquals(2.4, sMin)) step=0.1;
        }
        book(_h_cTheta,ihist,1,1);
      }
      else {
        throw Error("Invalid CMS energy for BELLE_2009_I822474");
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
        if (p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        if (ncount !=1 ) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (val.first==PID::PI0) {
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
        if (matched) {
          double cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
          if (cTheta<=0.8)    _sigmaEtaPi->fill(sqrtS());
          if (_h_cTheta ) _h_cTheta ->fill(cTheta);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_cTheta ) scale(_h_cTheta ,fact);
      scale(_sigmaEtaPi,fact);
      Estimate1DPtr tmp;
      book(tmp,1, 1, 1);
      barchart(_sigmaEtaPi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigmaEtaPi, _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I822474);

}
