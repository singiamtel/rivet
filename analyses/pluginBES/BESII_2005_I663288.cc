// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> phi pi+pi- K+K-
  class BESII_2005_I663288 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2005_I663288);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(333);
      declare(PSI, "PSI");
      // histograms
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h_mass[ix][iy], 1+ix, 1, 1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // J/psi decay products
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      if (PSI.decaying().size()!=1) vetoEvent;
      unsigned int imode=0;
      if      (PSI.modeMatches(0,3,mode1)) imode=0;
      else if (PSI.modeMatches(0,3,mode2)) imode=1;
      else  vetoEvent;
      unsigned int iMeson = imode == 0 ? 211 : 321;
      const Particle & Mp  = PSI.decayProducts()[0].at( iMeson)[0];
      const Particle & Mm  = PSI.decayProducts()[0].at(-iMeson)[0];
      const Particle & phi = PSI.decayProducts()[0].at(    333)[0];
      // mass histograms
      FourMomentum pX = Mp.momentum()+Mm.momentum();
      const double mX = pX.mass();
      _h_mass[imode][0]->fill(mX/GeV);
      _h_mass[imode][1]->fill((Mp.momentum()+phi.momentum()).mass()/GeV);
      _h_mass[imode][1]->fill((Mm.momentum()+phi.momentum()).mass()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0; ix<2; ++ix) {
        normalize(_h_mass[ix], 1.0, false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2][2];
    const map<PdgId,unsigned int> mode1 = { { 333,1}, { 211,1},{-211,1}};
    const map<PdgId,unsigned int> mode2 = { { 333,1}, { 321,1},{-321,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2005_I663288);

}
