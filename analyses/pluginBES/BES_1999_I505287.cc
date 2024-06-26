// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> K*0 Kbar*0
  class BES_1999_I505287 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1999_I505287);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::K0S);
      PSI.addStable(PID::PI0);
      PSI.addStable( 313);
      PSI.addStable(-313);
      declare(PSI, "PSI");
      // histograms
      book(_h_mass,1,1,1);
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_angle[ix],2,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for (unsigned int ix=0;ix<PSI.decaying().size();++ix) {
        if (!PSI.modeMatches(ix,3,mode)) continue;
        const Particle  & Kstar0  = PSI.decayProducts()[0].at( 313)[0];
        const Particle  & KstarB0 = PSI.decayProducts()[0].at(-313)[0];
        FourMomentum phad = Kstar0.mom()+KstarB0.mom();
        _h_mass->fill(phad.mass());
        Particle Kp;
        if (Kstar0.children()[0].pid()== 321 &&
            Kstar0.children()[1].pid()==-211) {
          Kp = Kstar0.children()[0];
        }
        else if (Kstar0.children()[1].pid()== 321 &&
                 Kstar0.children()[0].pid()==-211) {
          Kp = Kstar0.children()[1];
        }
        else {
          continue;
        }
        Particle Km;
        if (KstarB0.children()[0].pid()==-321 &&
            KstarB0.children()[1].pid()== 211) {
          Km = KstarB0.children()[0];
        }
        else if (KstarB0.children()[1].pid()==-321 &&
                 KstarB0.children()[0].pid()== 211) {
          Km = KstarB0.children()[1];
        }
        else {
          continue;
        }
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(phad.betaVec());
        // K*0 decay
        FourMomentum pKstar0 = boost.transform(Kstar0 .mom());
        Vector3 axis1 = pKstar0.p3().unit();
        const LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(pKstar0 .betaVec());
        FourMomentum pKp = boost1.transform(boost.transform(Kp.mom()));
        double cTheta1 = pKp.p3().unit().dot(axis1);
        _h_angle[0]->fill(cTheta1);
        Vector3 trans1 = pKp.p3() - cTheta1*pKp.p3().mod()*axis1;
        // Kbar*0 decay
        FourMomentum pKstarB0 = boost.transform(KstarB0.mom());
        Vector3 axis2 = pKstar0.p3().unit();
        const LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKstarB0.betaVec());
        FourMomentum pKm = boost2.transform(boost.transform(Km.mom()));
        double cTheta2 = pKm.p3().unit().dot(axis2);
        _h_angle[0]->fill(cTheta2);
        Vector3 trans2 = pKm.p3() - cTheta1*pKm.p3().mod()*axis2;
        double chi = abs(atan2(trans1.cross(trans2).dot(axis1),trans1.dot(trans2)))/M_PI*180.;
        _h_angle[1]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      normalize(_h_angle, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass, _h_angle[2];
    const map<PdgId,unsigned int> mode = { { 313,1}, {-313,1}, { 22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_1999_I505287);

}
