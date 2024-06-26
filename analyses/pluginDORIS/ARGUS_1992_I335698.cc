// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> pi- pi0 nu_tau
  class ARGUS_1992_I335698 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1992_I335698);


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
      // hists
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
      	if (!(TAU.modeMatches(ix,3,mode) || TAU.modeMatches(ix,3,modeCC))) continue;
        const Particle& pi0 = TAU.decayProducts()[ix].at( 111     )[0];
        const Particle& pim = TAU.decayProducts()[ix].at(-211*sign)[0];
        FourMomentum phad = pim.mom()+pi0.mom();
        _h[0]->fill(phad.mass());
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(phad.betaVec());
        FourMomentum ppi = boost.transform(pim.mom());
        const double cTheta = ppi.p3().unit().dot(phad.p3().unit());
        _h[1]->fill(cTheta);
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
    const map<PdgId,unsigned int> mode   = { { 111,1},{-211,1},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 111,1},{ 211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1992_I335698);

}
