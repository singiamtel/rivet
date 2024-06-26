// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> eta' rho+
  class CLEOII_1998_I466173 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1998_I466173);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::RHOPLUS);
      DS.addStable(PID::ETAPRIME);
      declare(DS,"DS");
      // histos
      book(_h, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for (unsigned int ix=0; ix < DS.decaying().size(); ++ix) {
        int sign = 1;
        if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,2,mode))         sign = 1;
        else if  (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,2,modeCC)) sign =-1;
        else continue;
        const Particle& rhoP = DS.decayProducts()[ix].at( sign*213)[0];
        // rho decay products
        Particle piP;
        if (rhoP.children().size()!=2) vetoEvent;
        if (rhoP.children()[0].pid()==211 &&
            rhoP.children()[1].pid()==111) {
          piP = rhoP.children()[0];
        }
        else if (rhoP.children()[1].pid()==211 &&
                 rhoP.children()[0].pid()==111) {
          piP = rhoP.children()[1];
        }
        else  vetoEvent;
        // compute the helicity angle
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(DS.decaying()[ix].mom().betaVec());
        FourMomentum prho=boost1.transform(rhoP.mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
        Vector3 axis = boost2.transform(boost1.transform(piP.mom())).p3().unit();
        _h->fill(axis.dot(prho.p3().unit()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    const map<PdgId,unsigned int> mode   = { { 213,1}, { 331,1}};
    const map<PdgId,unsigned int> modeCC = { {-213,1}, { 331,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1998_I466173);

}
