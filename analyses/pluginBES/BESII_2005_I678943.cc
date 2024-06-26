// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief $\pi^+\pi^-$ mass distribution in $J/\psi\to 2\pi^+2\pi^-$
  class BESII_2005_I678943 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2005_I678943);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      declare(PSI, "PSI");
      // histos
      book(_h, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      for (unsigned int ix=0; ix<PSI.decaying().size(); ++ix) {
        if (!PSI.modeMatches(ix,4,mode)) continue;
        const Particles& pip = PSI.decayProducts()[0].at( 211);
        const Particles& pim = PSI.decayProducts()[0].at(-211);
        bool kVeto=false;
        vector<double> mpipi; mpipi.reserve(4);
        for (unsigned int iy=0; iy<2; ++iy) {
          for (unsigned int iy=0; iy<2; ++iy) {
            double mpp = (pip[ix].mom()+pim[0].mom()).mass();
            kVeto |= (mpp>0.47&&mpp<0.53);
            mpipi.push_back(mpp);
          }
        }
        if (kVeto) continue;
        for (const double& m : mpipi) _h->fill(m);
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
    const map<PdgId,unsigned int> mode = { { 211,2}, { -211,2}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2005_I678943);

}
