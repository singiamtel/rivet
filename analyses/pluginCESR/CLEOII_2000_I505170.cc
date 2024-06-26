// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  tau -> 4 pions
  class CLEOII_2000_I505170 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2000_I505170);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      TAU.addStable(221);
      declare(TAU, "TAU");
      // histos
      for (unsigned int ix=0; ix<6; ++ix) {
        book(_h_mass[ix], 3, 1, 1+ix);
        if (ix>2) continue;
        book(_h_spect[ix], 1, 1, 1+ix);
      }
      book(_h_angle, 2, 1, 1);
      book(_c,"TMP/nTau");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
      	_c->fill();
        int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
        if (!(TAU.modeMatches(ix,5,mode  ) ||
              TAU.modeMatches(ix,5,modeCC))) continue;
        const Particle  & pi0 = TAU.decayProducts()[ix].at( 111     )[0];
        const Particles & pim = TAU.decayProducts()[ix].at(-211*sign);
        const Particle  & pip = TAU.decayProducts()[ix].at( 211*sign)[0];
        FourMomentum phad = pim[0].mom()+pim[1].mom()+pip.mom()+pi0.mom();
        const double q = phad.mass();
        const double mtau = TAU.decaying()[ix].mass();
        _h_spect[0]->fill(phad.mass(),sqr(mtau)*mtau/q/sqr(sqr(mtau)-sqr(q))/(sqr(mtau)+2.*sqr(q)));
        Particle omega;
        int iomega=-1;
        for (unsigned int ix=0; ix<2; ++ix) {
          Particle parent = pim[ix];
          while (!parent.parents().empty()) {
            parent = parent.parents()[0];
            if (parent.pid()==223) break;
          }
          if (parent.pid()==223) {
            omega=parent;
            iomega=ix;
          }
        }
        if (iomega>=0) {
          _h_spect[1]->fill(phad.mass(),sqr(mtau)*mtau/q/sqr(sqr(mtau)-sqr(q))/(sqr(mtau)+2.*sqr(q)));
          const int ibatch = iomega == 0 ? 1 : 0;
          const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(omega.mom().betaVec());
          const Vector3 axis = boost.transform(pim[ibatch].mom());
          const FourMomentum ppi0 = boost.transform(pi0.mom());
          const FourMomentum ppip = boost.transform(pip.mom());
          const Vector3 norm = ppi0.p3().cross(pip.p3()).unit();
          const double cTheta = axis.dot(norm);
          _h_angle->fill(cTheta);
        }
        else {
          _h_spect[2]->fill(phad.mass(),sqr(mtau)*mtau/q/sqr(sqr(mtau)-sqr(q))/(sqr(mtau)+2.*sqr(q)));
        }
        _h_mass[0]->fill((phad-pim[0].mom()).mass());
        _h_mass[0]->fill((phad-pim[1].mom()).mass());
        _h_mass[1]->fill((pim[0].mom()+pi0.mom()).mass());
        _h_mass[1]->fill((pim[1].mom()+pi0.mom()).mass());
        _h_mass[2]->fill((pim[1].mom()+pip.mom()).mass());
        _h_mass[3]->fill((pip   .mom()+pi0.mom()).mass());
        _h_mass[4]->fill(q);
        _h_mass[5]->fill((phad-pi0.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double Gamma = 6.582119514e-25/290.3e-15;
      const double GF    = 1.1663787e-5;
      const double cosC  = 0.97373;
      scale(_h_spect, Gamma*16.*sqr(M_PI)/sqr(GF*cosC)/ *_c);
      normalize(_h_mass,  1.0, false);
      normalize(_h_angle, 1.0, false);

    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect[3],_h_angle,_h_mass[6];
    CounterPtr _c;
    const map<PdgId,unsigned int> mode   = { { 111,1},{-211,2},{ 211,1},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 111,1},{ 211,2},{-211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2000_I505170);

}
