// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> K+K- mu+ nu_mu
  class BESIII_2023_I2674768 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2674768);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      DS.addStable(PID::ETA);
      DS.addStable(PID::ETAPRIME);
      declare(DS, "DS");
      // Book histograms
      for (unsigned int ix=0; ix<5; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for (unsigned int ix=0; ix<DS.decaying().size(); ++ix) {
        if (!DS.modeMatches(ix, 4, mode)) continue;
        const Particle& Kp = DS.decayProducts()[ix].at( 321)[0];
        const Particle& Km = DS.decayProducts()[ix].at(-321)[0];
        const Particle& ep = DS.decayProducts()[ix].at( -13)[0];
        const Particle& nue= DS.decayProducts()[ix].at(  14)[0];
        FourMomentum pPhi = Kp.mom()+Km.mom();
        _h[0]->fill(pPhi.mass());
        FourMomentum qq = DS.decaying()[ix].mom()-pPhi;
        _h[1]->fill(qq.mass2());
        // boost momenta to DS rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(DS.decaying()[ix].mom().betaVec());
        FourMomentum pPHI = boost.transform(pPhi);
        Matrix3 ptoz(-pPHI.p3().unit(), Vector3(0,0,1));
        boost.preMult(ptoz);
        // the momenta in frane to W along z
        FourMomentum pD  = boost.transform(DS.decaying()[ix].mom());
        FourMomentum pKp = boost.transform(Kp .mom());
        FourMomentum pKm = boost.transform(Km .mom());
        FourMomentum pe  = boost.transform(ep .mom());
        FourMomentum pnu = boost.transform(nue.mom());
        pPhi = pKp+pKm;
        qq = pD-pPhi;
        LorentzTransform boostK = LorentzTransform::mkFrameTransformFromBeta(pPhi.betaVec());
        Vector3 axisK = boostK.transform(pKp).p3().unit();
        _h[2]->fill(axisK.dot(pPhi.p3().unit()));
        LorentzTransform boostW = LorentzTransform::mkFrameTransformFromBeta(  qq.betaVec());
        Vector3 axisE = boostW.transform(pe).p3().unit();
        _h[3]->fill(axisE.dot(qq.p3().unit()));
        axisK.setZ(0.);
        axisE.setZ(0.);
        const double chi = atan2(axisE.cross(axisK).dot(qq.p3().unit()), axisE.dot(axisK));
        _h[4]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[5];
    const map<PdgId,unsigned int> mode = { { 321,1}, {-321,1}, {-13,1}, { 14,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2674768);

}
