// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> omega pi+ pi-
  class BESII_2004_I652399 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2004_I652399);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(223);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // histogram
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass[ix],1,1,1+ix);
        for (unsigned int iy=0; iy<4; ++iy) {
          book(_h_angle[iy][ix],2,1+iy,1+ix);
        }
      }
    }

    void findChildren(const Particle& p, Particles& pim, Particles& pip,
                      Particles& pi0, unsigned int &ncount) {
      for (const Particle& child : p.children()) {
        if (child.pid()==PID::PIPLUS) {
          pip.push_back(child);
          ncount+=1;
        }
        else if (child.pid()==PID::PIMINUS) {
          pim.push_back(child);
          ncount+=1;
        }
        else if (child.pid()==PID::PI0) {
          pi0.push_back(child);
          ncount+=1;
        }
        else if (child.children().empty()) {
          ncount+=1;
        }
        else {
          findChildren(child,pim,pip,pi0,ncount);
        }
      }
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
      // J/psi decay products
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if ( PSI.decaying().size()!=1) vetoEvent;
      if (!PSI.modeMatches(0,3,mode)) vetoEvent;
      const Particle& Pip   = PSI.decayProducts()[0].at( 211)[0];
      const Particle& Pim   = PSI.decayProducts()[0].at(-211)[0];
      const Particle& omega = PSI.decayProducts()[0].at( 223)[0];
      // mass histograms
      FourMomentum pPiPi = Pip.mom()+Pim   .mom();
      const double mPiPi = pPiPi.mass();
      _h_mass[0]->fill(mPiPi/GeV);
      _h_mass[1]->fill((Pip.mom()+omega.mom()).mass()/GeV);
      _h_mass[1]->fill((Pim.mom()+omega.mom()).mass()/GeV);
      unsigned int iplot=0;
      if      (mPiPi>.3 && mPiPi<.5) iplot=0;
      else if (mPiPi>.5 && mPiPi<.7) iplot=1;
      else  return;
      if (abs(Pip.p3().unit().dot(axis))>0.84) return;
      if (abs(Pim.p3().unit().dot(axis))>0.84) return;
      LorentzTransform boost;
      if (PSI.decaying()[0].mom().p3().mod() > 1*MeV) {
        boost = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].mom().betaVec());
      }
      pPiPi = boost.transform(pPiPi);
      // omega decay
      unsigned int ncount=0;
      Particles pip, pim, pi0;
      findChildren(omega,pim,pip,pi0,ncount);
      if ( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1))  return;
      if (abs(pip[0].p3().unit().dot(axis))>0.84)  return;
      if (abs(pim[0].p3().unit().dot(axis))>0.84)  return;
      FourMomentum pOmega = boost.transform(omega.mom());
      Vector3 e1Z = pOmega.p3().unit();
      Vector3 e1Y = e1Z.cross(axis).unit();
      Vector3 e1X = e1Y.cross(e1Z).unit();
      const double cOmega = e1Z.dot(axis);
      _h_angle[1][iplot]->fill(cOmega);
      // pion angles
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pPiPi.betaVec());
      FourMomentum pPi = boost2.transform(boost.transform(Pip.mom()));
      const double cPi = pPi.p3().unit().dot(e1Z);
      _h_angle[2][iplot]->fill(cPi);
      Vector3 trans1 = pPi.p3().unit()-cPi*e1Z;
      // omega angles
      LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
      FourMomentum ppip = boost3.transform(boost.transform(pip[0].mom()));
      FourMomentum ppim = boost3.transform(boost.transform(pim[0].mom()));
      Vector3 nW = ppip.p3().cross(ppim.p3()).unit();
      const double bW = nW.dot(e1X);
      _h_angle[3][iplot]->fill(bW);
      Vector3 trans2 = ppip.p3()-ppip.p3().dot(e1Z)*e1Z;
      const double chi = abs(atan2(trans1.cross(trans2).dot(e1Z),trans1.dot(trans2)));
      _h_angle[0][iplot]->fill(chi/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      for (unsigned int iy=0; iy<4; ++iy) {
        normalize(_h_angle[iy], 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_angle[4][2];
    const map<PdgId,unsigned int> mode = { { 223,1}, { 211,1},{-211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2004_I652399);

}
