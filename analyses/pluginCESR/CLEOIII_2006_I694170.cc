// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Upsilon -> gamma pi+pi- K+K- ppbar
  class CLEOIII_2006_I694170 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOIII_2006_I694170);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==553);
      declare(ufs, "UFS");
      DecayedParticles UPS(ufs);
      declare(UPS, "UPS");
      declare(Beam(), "Beams");
      // histos
      for (unsigned int ix=0;ix<4;++ix) {
        book(_h_mass[ix],1+ix,1,1);
      }
      for (unsigned int ix=0;ix<2;++ix) {
        for (unsigned int iy=0;iy<2;++iy) {
	        book(_h_angle[ix][iy],5+ix,1,1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.pid()>0) {
        axis = beams.first .momentum().p3().unit();
      }
      else {
        axis = beams.second.momentum().p3().unit();
      }
      unsigned int iprod[3]={211,321,2212};
      const map<PdgId,unsigned int> modes[3] = {{ { 211 ,1}, { -211 ,1}, { 22,1}},
						{ { 321 ,1}, { -321 ,1}, { 22,1}},
						{ { 2212,1}, { -2212,1}, { 22,1}}};
      DecayedParticles UPS = apply<DecayedParticles>(event, "UPS");
      if (UPS.decaying().size()!=1) vetoEvent;
      int imode=-1;
      for (unsigned int ix=0; ix<3; ++ix) {
	      if (UPS.modeMatches(0,3,modes[ix])) {
          imode=ix;
          break;
        }
      }
      if (imode<0) vetoEvent;
      const Particle & Mp  = UPS.decayProducts()[0].at( iprod[imode])[0];
      const Particle & Mm  = UPS.decayProducts()[0].at(-iprod[imode])[0];
      FourMomentum pMM = Mp.momentum()+Mm.momentum();
      double mass = pMM.mass();
      if (imode==0) {
        _h_mass[0]->fill(mass);
        _h_mass[1]->fill(mass);
      }
      else {
        _h_mass[1+imode]->fill(mass);
      }
      if (imode==2) return;
      bool hasf2[2]={false,false};
      for (const Particle & p : UPS.decaying()[0].children()) {
        if (p.pid()==225)      hasf2[0]=true;
        else if (p.pid()==335) hasf2[1]=true;
      }
      if ( (imode==0 && hasf2[0]) || (imode==1&&hasf2[1]) ) {
        const Particle & gam = UPS.decayProducts()[0].at(22)[0];
        const double cTheta = axis.dot(gam.p3().unit());
        _h_angle[imode][1]->fill(abs(cTheta));
        // remaining angles
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(UPS.decaying()[0].mom().betaVec());
        FourMomentum pGamma = boost1.transform(gam.mom());
        pMM = boost1.transform(pMM);
        Vector3 e1z = pGamma.p3().unit();
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pMM.betaVec());
        Vector3 axis2 = boost2.transform(boost1.transform(Mp.mom())).p3().unit();
        _h_angle[imode][0]->fill(abs(e1z.dot(axis2)));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<4; ++ix) {
        normalize(_h_mass[ix], 1.0, false);
      }
      for (unsigned int ix=0; ix<2; ++ix) {
	      for (unsigned int iy=0; iy<2; ++iy) {
	        normalize(_h_angle[ix][iy], 1.0, false);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[4], _h_angle[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOIII_2006_I694170);

}
