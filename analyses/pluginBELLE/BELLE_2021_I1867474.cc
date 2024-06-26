// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> pi+pi-eta, K+K-eta
  class BELLE_2021_I1867474 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1867474);


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
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<3; ++iy) {
          book(_h[ix][iy],1+ix,1,1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for (unsigned int ix=0; ix<D0.decaying().size(); ++ix) {
        unsigned int imode=0;
        if       (D0.modeMatches(ix,3,mode1)) imode=0;
        else if  (D0.modeMatches(ix,3,mode2)) imode=1;
        else continue;
        int sign = D0.decaying()[ix].pid()/D0.decaying()[ix].abspid();
        int im = 211+imode*110;
        const Particle& eta = D0.decayProducts()[ix].at(221)[0];
        const Particle& pim = D0.decayProducts()[ix].at(-sign*im)[0];
        const Particle& pip = D0.decayProducts()[ix].at( sign*im)[0];
        _h[imode][0]->fill((pip.mom()+pim.mom()).mass2());
        _h[imode][1]->fill((pim.mom()+eta.mom()).mass2());
        _h[imode][2]->fill((pip.mom()+eta.mom()).mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h[ix], 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    const map<PdgId,unsigned int> mode1 = { { 211,1},{-211,1}, {221,1}};
    const map<PdgId,unsigned int> mode2 = { { 321,1},{-321,1}, {221,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1867474);

}
