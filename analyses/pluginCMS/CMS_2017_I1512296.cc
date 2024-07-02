// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi at 5.02 TeV
  class CMS_2017_I1512296 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1512296);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // histograms
      for (int ix=0; ix<2; ++ix) {
        book(_h_JPsi_pT[ix], {-2.4,-1.93,-1.5,-0.9,0.,0.9,1.5,1.93,2.4});
        for (int iy=0; iy<8; ++iy) {
          Histo1DPtr tmp;
          if (iy<4) {
            book(_h_JPsi_pT[ix]->bin(iy+1), 2+13*ix, 1, 4-iy);
          }
          else {
            book(_h_JPsi_pT[ix]->bin(iy+1), 1+13*ix, 1, 8-iy);
          }
        }
      }
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_JPsi_y[ix], {6.5,10.,30.});
        for (auto& b : _h_JPsi_y[ix]->bins()) {
          book(b, 5+13*ix, 1, b.index());
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
        // prompt/non-prompt
        bool nonPrompt = p.fromBottom();
        const double rap = p.rap();
        const double xp = p.perp();
        _h_JPsi_pT[nonPrompt]->fill(rap,xp);
        _h_JPsi_y [nonPrompt]->fill(xp,rap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br to muons PDG 2021
      const double br = 0.05961;
      const double factor = br*crossSection()/microbarn/sumOfWeights();
      const vector<double> pTbins={6.5,10.,30.};
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_JPsi_pT[ix], factor);
        divByGroupWidth(_h_JPsi_pT[ix]);
        // in y just single differential, don't divide by group width
        scale(_h_JPsi_y[ix], factor);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_JPsi_pT[2],_h_JPsi_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1512296);

}
