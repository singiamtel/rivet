// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> rho pi
  class BELLE_2002_I589587 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2002_I589587);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511 || Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable( 111);
      declare(BB, "BB");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass [ix], 1, 1, 1+ix);
        book(_h_angle[ix], 2, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0; ix<BB.decaying().size(); ++ix) {
      	int sign = BB.decaying()[ix].pid()/BB.decaying()[ix].abspid();
        unsigned int imode=0;
      	if (BB.modeMatches(ix,3,mode1)) imode=0;
        else if(sign== 1 && BB.modeMatches(ix,3,mode2  )) imode=1;
        else if(sign==-1 && BB.modeMatches(ix,3,mode2CC)) imode=1;
        else continue;
        const Particles & pip = BB.decayProducts()[ix].at( sign*211);
        const Particles & pim = BB.decayProducts()[ix].at(-sign*211);
        bool pipiVeto=false;
        Particle pRho[2];
        double delta=1e30;
        for (const Particle& p1 : pip) {
          for(const Particle& p2 : pim) {
            double mpipi = (p1.mom()+p2.mom()).mass();
            if (abs(mpipi-1.86484)<0.14 ||
                abs(mpipi-3.09690)<0.07 ||
                abs(mpipi-3.68610)<0.05) {
              pipiVeto=true;
              break;
            }
            if (abs(mpipi-.77)<delta) {
              delta= abs(mpipi-.77);
              pRho[0]=p1;
              pRho[1]=p2;
            }
          }
        }
        if (pipiVeto) continue;
        if (imode==0) {
          const Particle & pi0 = BB.decayProducts()[ix].at(111)[0];
          delta=1e30;
          for (const Particle& p1 : pip) {
            double mpipi = (p1.mom()+pi0.mom()).mass();
            if (abs(mpipi-1.86966)<0.05) {
              pipiVeto=true;
              break;
            }
            if (abs(mpipi-.77)<delta) {
              delta= abs(mpipi-.77);
              pRho[0]=p1;
              pRho[1]=pi0;
            }
          }
          if (pipiVeto) continue;
          for (const Particle& p1 : pim) {
            double mpipi = (p1.mom()+pi0.mom()).mass();
            if (abs(mpipi-1.86966)<0.05) {
              pipiVeto=true;
              break;
            }
            if (abs(mpipi-.77)<delta) {
              delta= abs(mpipi-.77);
              pRho[0]=p1;
              pRho[1]=pi0;
            }
          }
          if (pipiVeto) continue;
        }
        FourMomentum prho = pRho[0].mom()+pRho[1].mom();
        double mpipi = prho.mass();
        _h_mass[imode]->fill(mpipi);
        // check if we should compute helicty angle
        if      (imode==0 && (mpipi<.6 || mpipi>0.95)) continue;
        else if (imode==1 && (mpipi<.62|| mpipi>0.92)) continue;
        // compute the helicity angle
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].mom().betaVec());
        prho = boost1.transform(prho);
        FourMomentum ppi  = boost1.transform(pRho[0].mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
        ppi = boost2.transform(ppi);
        double cTheta = ppi.p3().unit().dot(prho.p3().unit());
        _h_angle[imode]->fill(cTheta);
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
    Histo1DPtr _h_mass[2], _h_angle[2];
    const map<PdgId,unsigned int> mode1   = { { 211,1}, {-211,1}, { 111,1} };
    const map<PdgId,unsigned int> mode2   = { { 211,2}, {-211,1} };
    const map<PdgId,unsigned int> mode2CC = { {-211,2}, { 211,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2002_I589587);

}
