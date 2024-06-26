// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Transversality angles for $\eta_c\to f_2f^{(\prime)}_2$
  class BELLE_2008_I754089 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2008_I754089);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==441);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(225);
      ETA.addStable(335);
      declare(ETA, "ETA");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for (unsigned int ix=0; ix<ETA.decaying().size(); ++ix) {
        int imode=-1;
      	if      (ETA.modeMatches(ix,2,mode1)) imode=0;
      	else if (ETA.modeMatches(ix,2,mode2)) imode=1;
      	else {
          continue;
        }
        Vector3 trans[2];
        bool found=true;
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(ETA.decaying()[ix].mom().betaVec());
        Vector3 axis = boost1.transform(ETA.decaying()[ix].children()[0].mom()).p3().unit();
        for (unsigned int iy=0;iy<ETA.decaying()[ix].children().size();++iy) {
          int iMeson = ETA.decaying()[ix].children()[iy].pid()==225 ? 211 : 321;
          if (ETA.decaying()[ix].children()[iy].children().size()!=2) found=false;
          if (ETA.decaying()[ix].children()[iy].children()[0].pid()!=
             -ETA.decaying()[ix].children()[iy].children()[1].pid()) found=false;
          if (ETA.decaying()[ix].children()[iy].children()[0].abspid()!=iMeson) found=false;
          if (!found) break;
          FourMomentum pF2 = boost1.transform(ETA.decaying()[ix].children()[iy].mom());
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pF2.betaVec());
          FourMomentum pMeson = boost2.transform(boost1.transform(ETA.decaying()[ix].children()[iy].children()[0].mom()));
          trans[iy] = pMeson.p3().unit()-pMeson.p3().unit().dot(axis)*axis;
        }
        if (!found) continue;
        const double chi = atan2(trans[0].cross(trans[1]).dot(axis),trans[0].dot(trans[1]));
        _h[imode]->fill(chi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode1 = { { 225,2}};
    const map<PdgId,unsigned int> mode2 = { { 225,1},{ 335,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2008_I754089);

}
