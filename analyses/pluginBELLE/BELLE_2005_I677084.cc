// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B- -> D0 D*-
  class BELLE_2005_I677084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I677084);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==521), "UFS");

      // Book histograms
      book(_h_cTheta,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over B0 mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        // find the B decay
        int sign = p.pid()/p.abspid();
        if (p.children().size()!=2) continue;
        Particle Dstar;
        if (p.children()[0].pid()== sign*413 &&
            p.children()[1].pid()==-sign*421) {
          Dstar = p.children()[0];
        }
        else if (p.children()[1].pid()== sign*413 &&
          p.children()[0].pid()==-sign*421) {
          Dstar = p.children()[1];
        }
        else {
          continue;
        }
        // and the D* decay
        if (Dstar.children().size()!=2) continue;
        Particle pip;
        if (Dstar.children()[0].pid()==sign*421 &&
            Dstar.children()[1].pid()==sign*211) {
          pip = Dstar.children()[1];
        }
        else if (Dstar.children()[1].pid()==sign*421 &&
                 Dstar.children()[0].pid()==sign*211) {
          pip = Dstar.children()[0];
        }
        else {
          continue;
        }
        // compute the helicity angle
        // boost to B+ rest frame
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        FourMomentum pDstar = boost1.transform(Dstar.momentum());
        FourMomentum pPi    = boost1.transform(pip  .momentum());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
        pPi = boost2.transform(pPi);
        const double cTheta = pPi.p3().unit().dot(pDstar.p3().unit());
        _h_cTheta->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cTheta);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I677084);

}
