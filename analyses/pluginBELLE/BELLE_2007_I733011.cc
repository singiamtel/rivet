// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Bbar0 -> D0 pi+ pi-
  class BELLE_2007_I733011 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I733011);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 421);
      B0.addStable(-421);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<3;++ix) {
      	book(_h_pipi[ix], 3,1,1+ix);
      	if (ix==2) continue;
      	book(_h_mass[ix], 1,1,1+ix);
      	book(_h_Dpi [ix], 2,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode))        sign= 1;
      	else if (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) sign=-1;
        else {
          continue;
        }
        const Particle& pip = B0.decayProducts()[ix].at( sign*211)[0];
        const Particle& pim = B0.decayProducts()[ix].at(-sign*211)[0];
        const Particle& D0  = B0.decayProducts()[ix].at(-sign*421)[0];
        // compute the helicity angles
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].mom().betaVec());
        // D pi mass and angle
        FourMomentum pDpi = D0.mom()+pim.mom();
        double mDpi = pDpi.mass();
        pDpi = boost1.transform(pDpi);
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDpi.betaVec());
        FourMomentum pD = boost2.transform(boost1.transform(D0.mom()));
        double cDpi = pD.p3().unit().dot(pDpi.p3().unit());
        // fill histos
        if      (abs(mDpi-2.46)<.1) _h_Dpi[0]->fill(cDpi);
        else if (abs(mDpi-2.30)<.1) _h_Dpi[1]->fill(cDpi);
        if (cDpi>0.) _h_mass[0]->fill(mDpi);
        // pi pi mass and angle
        FourMomentum ppipi = pim.mom()+pip.mom();
        double mpipi = ppipi.mass();
        ppipi = boost1.transform(ppipi);
        LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
        FourMomentum ppip = boost3.transform(boost1.transform(pip.mom()));
        double cpipi = -ppip.p3().unit().dot(ppipi.p3().unit());
        // fill angles
        if      (abs(mpipi-0.78)<0.2) _h_pipi[0]->fill(cpipi);
        else if (abs(mpipi-1.20)<0.1) _h_pipi[1]->fill(cpipi);
        else if (mpipi<0.6)           _h_pipi[2]->fill(cpipi);
        if (cpipi>0.) _h_mass[1]->fill(mpipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi, 1.0, false);
      normalize(_h_mass, 1.0, false);
      normalize(_h_Dpi,  1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_Dpi[2],_h_pipi[3];
    const map<PdgId,unsigned int> mode   = { { 211,1}, {-211,1}, {-421,1}};
    const map<PdgId,unsigned int> modeCC = { { 211,1}, {-211,1}, { 421,1}};
    /// @}


      };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I733011);

}
