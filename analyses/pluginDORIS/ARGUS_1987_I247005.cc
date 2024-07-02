// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 phi
  class ARGUS_1987_I247005 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I247005);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::PHI);
      declare(D0, "D0");
      // histos
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for (unsigned int ix=0; ix<D0.decaying().size(); ++ix) {
        if (!D0.modeMatches(ix,2,mode)) continue;
        //const Particle & K0  = D0.decayProducts()[ix].at(310)[0];
        const Particle& phi = D0.decayProducts()[ix].at(333)[0];
        if (phi.children().size()!=2) continue;
        if (phi.children()[0].abspid()!=PID::KPLUS ||
            phi.children()[1].abspid()!=PID::KPLUS ) continue;
        // first boost all relevant momenta to D0 rest frame
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(D0.decaying()[ix].mom());
        FourMomentum pphi = boost.transform(phi.mom());
        FourMomentum pK   = boost.transform(phi.children()[0].mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pphi.betaVec());
        Vector3 axis1 = pphi.p3().unit();
        Vector3 axis2 = boost2.transform(pK ).p3().unit();
        _h->fill(axis1.dot(axis2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    const map<PdgId,unsigned int> mode = { { 310,1}, { 333,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I247005);

}
