// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi at 7 TeV
  class CMS_2011_I878118 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I878118);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h_psi[ix], {0.,1.2,1.6,2.4});
        for (unsigned int iy=1; iy<4; ++iy) {
          if (ix<2) {
            book(_h_psi[ix]->bin(iy), 1+iy+ix*9, 1, 1);
          }
          else {
            book(_h_psi[ix]->bin(iy), "TMP/total"+toString(iy), refData(4+iy,1,1));
          }
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for (const Particle & p : ufs.particles(Cuts::pid==443)) {
        // cuts on rapidity
        const double y = p.absrap();
        if (y>2.4) continue;
        const double pT = p.perp();
        // prompt
        unsigned int iprompt = p.fromBottom();
        _h_psi[iprompt]->fill(y,pT);
        _h_psi[   2   ]->fill(y,pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // branching ratio
      const double br = 0.05961;
      // 0.5 due folded rapidity
      scale(_h_psi, 0.5*br*crossSection() / nanobarn/ sumOfWeights());
      divByGroupWidth(_h_psi);
      // non-prompt fraction
      for (unsigned int iy=0; iy<_h_psi[1]->numBins(); ++iy) {
        Estimate1DPtr tmp;
        book(tmp, 5+iy, 1, 1);
        efficiency(_h_psi[1]->bin(iy+1), _h_psi[2]->bin(iy+1),tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_psi[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2011_I878118);

}
