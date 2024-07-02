// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> omega K+K-
  class BESII_2004_I658084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2004_I658084);


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
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_mass[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<4;++ix)
	book(_h_angle[ix],2,1,1+ix);
    }

    void findChildren(const Particle& p, Particles& pim, Particles& pip,
                      Particles & pi0, unsigned int &ncount) {
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
      const Particle& Kp    = PSI.decayProducts()[0].at( 321)[0];
      const Particle& Km    = PSI.decayProducts()[0].at(-321)[0];
      const Particle& omega = PSI.decayProducts()[0].at( 223)[0];
      // mass histograms
      FourMomentum pKK = Kp.mom()+Km   .mom();
      double mKK = pKK.mass();
      _h_mass[0]->fill(mKK/GeV);
      _h_mass[1]->fill((Kp.mom()+omega.mom()).mass()/GeV);
      _h_mass[1]->fill((Km.mom()+omega.mom()).mass()/GeV);
      if (mKK<2.) return;
      if (abs(Kp.p3().unit().dot(axis))>0.84) return;
      if (abs(Km.p3().unit().dot(axis))>0.84) return;
      LorentzTransform boost;
      if (PSI.decaying()[0].mom().p3().mod() > 1*MeV) {
        boost = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].mom().betaVec());
      }
      pKK = boost.transform(pKK);
      // omega decay
      unsigned int ncount=0;
      Particles pip,pim,pi0;
      findChildren(omega,pim,pip,pi0,ncount);
      if ( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) return;
      if (abs(pip[0].p3().unit().dot(axis))>0.84) return;
      if (abs(pim[0].p3().unit().dot(axis))>0.84) return;
      FourMomentum pOmega = boost.transform(omega.mom());
      Vector3 e1Z = pOmega.p3().unit();
      Vector3 e1Y = e1Z.cross(axis).unit();
      Vector3 e1X = e1Y.cross(e1Z).unit();
      const double cOmega = e1Z.dot(axis);
      _h_angle[1]->fill(cOmega);
      // kaon angles
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
      FourMomentum pK = boost2.transform(boost.transform(Kp.mom()));
      const double cK = pK.p3().unit().dot(e1Z);
      _h_angle[2]->fill(cK);
      Vector3 transK = pK.p3().unit()-cK*e1Z;
      // omega angles
      LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
      FourMomentum ppip = boost3.transform(boost.transform(pip[0].mom()));
      FourMomentum ppim = boost3.transform(boost.transform(pim[0].mom()));
      Vector3 nW = ppip.p3().cross(ppim.p3()).unit();
      const double bW = nW.dot(e1X);
      _h_angle[3]->fill(bW);
      Vector3 transPi = ppip.p3()-ppip.p3().dot(e1Z)*e1Z;
      const double chi = abs(atan2(transK.cross(transPi).dot(e1Z),transK.dot(transPi)));
      _h_angle[0]->fill(chi/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      normalize(_h_angle, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_angle[4];
    const map<PdgId,unsigned int> mode = { { 223,1}, { 321,1},{-321,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2004_I658084);

}
