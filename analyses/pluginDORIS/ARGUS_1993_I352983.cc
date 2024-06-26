// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+ pi-
  class ARGUS_1993_I352983 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_I352983);


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
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for (unsigned int ix=0; ix<D0.decaying().size(); ++ix) {
        if (!D0.modeMatches(ix,3,mode) ) continue;
        int sign = D0.decaying()[ix].pid()/421;
        const Particles & pip= D0.decayProducts()[ix].at( sign*211);
        const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
        const Particles & K0 = D0.decayProducts()[ix].at( 310);
        double mminus = (pim[0].mom()+K0[0].mom() ).mass2();
        double mplus  = (pip[0].mom()+K0[0].mom() ).mass2();
        double mpipi  = (pip[0].mom()+pim[0].mom()).mass2();
        _h[2]->fill(mplus);
        _h[1]->fill(mminus);
        _h[0]->fill(mpipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode = { { 310,1}, { 211,1},{-211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1993_I352983);

}
