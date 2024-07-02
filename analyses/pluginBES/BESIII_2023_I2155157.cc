// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma KS0 KS0 pi0
  class BESIII_2023_I2155157 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2155157);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::K0S);
      PSI.addStable(PID::ETA);
      PSI.addStable(PID::PI0);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // Book histograms
      for (unsigned int ix=0; ix<7; ++ix) {
      	book(_h[ix], 1, 1, 1+ix);
      }
    }

    // angle cuts due regions of BES calorimeter
    bool vetoPhoton(const double & cTheta) {
      return cTheta>0.92 || (cTheta>0.8 && cTheta<0.86);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
      	axis = beams.first .mom().p3().unit();
      }
      else {
      	axis = beams.second.mom().p3().unit();
      }
      // find the J/psi decays
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if ( PSI.decaying().size()!=1) vetoEvent;
      if (!PSI.modeMatches(0,4,mode)) vetoEvent;
      const Particle & pi0 = PSI.decayProducts()[0].at(111)[0];
      const Particles& K0  = PSI.decayProducts()[0].at(310);
      const Particle & gam = PSI.decayProducts()[0].at( 22)[0];
      double cTheta = axis.dot(gam.p3().unit());
      // photon angle cut
      if (vetoPhoton(abs(cTheta))) vetoEvent;
      // gamma pi omega mass cut
      if (abs((pi0.mom()+gam.mom()).mass()-0.78266)<0.04) vetoEvent;
      double mKKpi = (K0[0].mom()+K0[1].mom()+pi0.mom()).mass();
      if (mKKpi>1.6) vetoEvent;
      _h[0]->fill(mKKpi);
      _h[1]->fill((K0[0].mom()+K0[1].mom()).mass());
      _h[1]->fill((K0[0].mom()+pi0  .mom()).mass());
      _h[2]->fill((K0[1].mom()+pi0  .mom()).mass());
      _h[3]->fill(cTheta);
      // remaining angles
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].mom().betaVec());
      FourMomentum pGamma    = boost1.transform(gam.mom());
      FourMomentum pHadron = boost1.transform(K0[0].mom()+K0[1].mom()+pi0.mom());
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pHadron.betaVec());
      Vector3 axis1 = pGamma.p3().unit();
      Vector3 axis2 = boost2.transform(boost1.transform(pi0.mom())).p3().unit();
      for (unsigned int ix=0; ix<2; ++ix) {
        axis2 = boost2.transform(boost1.transform(K0[ix].mom())).p3().unit();
        _h[4]->fill(axis1.dot(axis2));
      }
      _h[5]->fill(axis1.dot(axis2));
      FourMomentum pKK = boost2.transform(boost1.transform(K0[0].mom()+K0[1].mom()));
      axis2=pKK.p3().unit();
      LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
      for (unsigned ix=0; ix<2; ++ix) {
        Vector3 axis3 = boost3.transform(boost2.transform(boost1.transform(K0[ix].mom()))).p3().unit();
        _h[6]->fill(axis3.dot(axis2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
     	normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[7];
    const map<PdgId,unsigned int> mode = { { 111,1}, { 310,2},{ 22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2155157);

}
