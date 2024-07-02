// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma eta pi+ pi-
  class BES_1999_I498114 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1999_I498114);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::ETA);
      declare(PSI, "PSI");
      // histos
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode = { { 221,1}, { 211,1}, { 211,1}, { 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for (unsigned int ix=0; ix<PSI.decaying().size(); ++ix) {
        if (!PSI.modeMatches(ix,4,mode)) continue;
        const Particle & gam = PSI.decayProducts()[0].at( 22)[0];
        _h->fill((PSI.decaying()[ix].mom()-gam.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_1999_I498114);

}
