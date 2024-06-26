// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> D0 Dbar0 K+
  class BELLE_2008_I756554 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I756554);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<int>("PID", 100433);
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 421);
      BP.addStable(-421);
      declare(BP, "BP");
      // histos
      for (unsigned int ix=0;ix<3;++ix) {
        book(_h_mass[ix],1,1,1+ix);
      }
      book(_h_angle,2,1,1);
      book(_c,"TMP/c");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
        int sign=1;
        if      (BP.modeMatches(ix,3,mode  )) sign= 1;
        else if (BP.modeMatches(ix,3,modeCC)) sign=-1;
        else continue;
        const Particle & Kp    = BP.decayProducts()[ix].at( sign*321)[0];
        const Particle & D0    = BP.decayProducts()[ix].at( sign*421)[0];
        const Particle & Dbar0 = BP.decayProducts()[ix].at(-sign*421)[0];
        _c->fill();
        FourMomentum pDD = D0.mom()+Dbar0.mom();
        double mDD = pDD.mass();
        _h_mass[0]->fill(mDD);
        FourMomentum pDK = D0.mom()+Kp.mom();
        if(mDD>3.85) _h_mass[2]->fill(pDK.mass());
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());
        pDD = boost1.transform(pDD);
        // boost to DD frame
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDD.betaVec());
        FourMomentum pD = boost2.transform(boost1.transform(D0.mom()));
        // helicity angle
        double cTheta = pDD.p3().unit().dot(pD.p3().unit());
        if (cTheta>0) _h_mass[1]->fill(mDD);
        // check for intermediate DS
        bool Ds2700 = false;
        for(const Particle & child : BP.decaying()[ix].children()) {
          Ds2700 |= child.abspid()==_pid;
        }
        if (!Ds2700) continue;
        pDK = boost1.transform(pDK);
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pDK.betaVec());
        FourMomentum pK = boost3.transform(boost1.transform(Kp.mom()));
        cTheta = pDK.p3().unit().dot(pK.p3().unit());
        _h_angle->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_mass, 1.0/ *_c);
      normalize(_h_angle);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_angle;
    CounterPtr _c;
    int _pid;
    const map<PdgId,unsigned int> mode   = { { 421,1}, {-421,1}, { 321,1}};
    const map<PdgId,unsigned int> modeCC = { { 421,1}, {-421,1}, {-321,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I756554);

}
