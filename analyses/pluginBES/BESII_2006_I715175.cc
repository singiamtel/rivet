// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma omega omega
  class BESII_2006_I715175 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I715175);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::OMEGA);
      declare(PSI, "PSI");
      declare(Beam(), "Beams");
      // histos
      for (unsigned int ix=0; ix<6; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }

    // angle cuts due regions of BES calorimeter
    bool vetoPhoton(const double& cTheta) {
      return cTheta>0.92 || (cTheta>0.8 && cTheta<0.86);
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
      Vector3 trans[2];
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
      if ( PSI.decaying().size()!=1)  vetoEvent;
      if (!PSI.modeMatches(0,3,mode)) vetoEvent;
      // particles
      const Particles& omega = PSI.decayProducts()[0].at(223);
      const Particle & gam   = PSI.decayProducts()[0].at( 22)[0];
      // photon polar angle
      double cTheta = axis.dot(gam.p3().unit());
      if (vetoPhoton(abs(cTheta))) vetoEvent;
      // remaining angles
      LorentzTransform boost1  = LorentzTransform::mkFrameTransformFromBeta(PSI.decaying()[0].mom().betaVec());
      FourMomentum pGamma      = boost1.transform(gam.mom());
      FourMomentum pOmegaOmega = boost1.transform(omega[0].mom()+omega[1].mom());
      Vector3 e1z = pGamma.p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pOmegaOmega.betaVec());
      Vector3 axis2 = boost2.transform(boost1.transform(omega[0].mom())).p3().unit();
      double cOmega[2], phiOmega[2], cN[2];
      for (unsigned int ix=0; ix<2; ++ix) {
        Vector3 axis3 = boost2.transform(boost1.transform(omega[ix].mom())).p3().unit();
        cOmega[ix] = e1z.dot(axis3) ;
        phiOmega[ix] = atan2(axis3.dot(e1y),axis3.dot(e1x));
        if (phiOmega[ix]<0.) phiOmega[ix]+=2.*M_PI;
        // omega decay
        unsigned int ncount=0;
        Particles pip,pim,pi0;
        findChildren(omega[ix],pim,pip,pi0,ncount);
        if (ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) vetoEvent;
        FourMomentum ppip = boost2.transform(boost1.transform(pip[0].mom()));
        FourMomentum ppim = boost2.transform(boost1.transform(pim[0].mom()));
        FourMomentum pOmega = boost2.transform(boost1.transform(omega[ix].mom()));
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
        // boost to omega frame
        Vector3 axisZ = pOmega.p3().unit();
        ppip = boost3.transform(ppip);
        ppim = boost3.transform(ppim);
        Vector3 norm = ppip.p3().cross(ppim.p3()).unit();
        cN[ix] = norm.dot(axisZ);
        trans[ix] = norm - cN[ix]*axisZ;
      }
      // finally fill all the histos
      // photon polar angle
      _h[0]->fill(cTheta);
      // azimuthal angle
      _h[1]->fill(gam.p3().phi()/M_PI*180.);
      // omega angles
      for (unsigned int ix=0; ix<2; ++ix) {
        _h[2]->fill(cOmega[ix]);
        _h[3]->fill(180.*phiOmega[ix]/M_PI);
        _h[4]->fill(cN[ix]);
      }
      const double chi = atan(trans[0].cross(trans[1]).dot(axis2)/trans[0].dot(trans[1]));
      _h[5]->fill(abs(chi)/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    const map<PdgId,unsigned int> mode = { { 223,2},{ 22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I715175);

}
