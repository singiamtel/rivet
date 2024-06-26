// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma pi pi
  class BESII_2006_I712991 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I712991);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::ETA);
      PSI.addStable(PID::PI0);
      declare(PSI, "PSI");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for (unsigned int ix=0; ix<PSI.decaying().size(); ++ix) {
        unsigned int imode=0;
        if (PSI.modeMatches(ix,3,mode1))      imode=0;
        else if (PSI.modeMatches(ix,3,mode2)) imode=1;
        else continue;
        const Particle& gam = PSI.decayProducts()[0].at( 22)[0];
        _h[imode]->fill((PSI.decaying()[ix].mom()-gam.mom()).mass());
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
    const map<PdgId,unsigned int> mode1= { { 211,1}, {-211,1}, { 22,1}};
    const map<PdgId,unsigned int> mode2= { { 111,2},{ 22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I712991);

}
