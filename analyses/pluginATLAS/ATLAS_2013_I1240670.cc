// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {


  /// @brief B+ meson at 7 TeV
  class ATLAS_2013_I1240670 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1240670);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(Cuts::abspid==521), "UFS");
      const vector<double> yedges = {0.0,0.5,1.0,1.5,2.25};
      book(_h_pT_y, yedges);
      for (auto& b : _h_pT_y->bins()) {
        book(b, b.index(),1,1);
      }
      book(_h_pT,5,1,1);
      book(_h_y ,6,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over B+ states
      for (const Particle& p : ufs.particles()) {
        // cuts on pT and rapidity
        double y  = p.absrap();
        if (y>2.25) continue;
        double pT = p.perp();
        _h_pT_y->fill(y,pT);
        _h_pT->fill(pT);
        if(pT>9. && pT<120.) _h_y->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br for B+-> J/psi K+ psi->mu+mu- (PDG2020)
      double br = 1.020e-3*0.05961;
      double fact = 0.5*br*crossSection()/picobarn/sumOfWeights();
      scale(_h_pT, fact);
      // 0.5 from +/- y
      scale(_h_pT_y, 0.5*fact);
      divByGroupWidth(_h_pT_y);
      // 0.5 from +/- y
      scale(_h_y, 0.5*fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_pT_y;
    Histo1DPtr _h_pT,_h_y;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1240670);

}
