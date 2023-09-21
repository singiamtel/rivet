// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// LHCb prompt charm hadron pT and rapidity spectra
  class LHCB_2013_I1218996 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1218996);

    /// @}


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      declare(UnstableParticles(), "UFS");

      /// Book histograms
      book(_h_pdg421_Dzero_pT_y,     {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg411_Dplus_pT_y,     {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg413_Dstarplus_pT_y, {2., 2.5, 3., 3.5, 4., 4.5});
      book(_h_pdg431_Dsplus_pT_y,    {2., 2.5, 3., 3.5, 4., 4.5});
      for (size_t i = 1; i < _h_pdg421_Dzero_pT_y->numBins()+1; ++i) {
        size_t y = _h_pdg421_Dzero_pT_y->bin(i).index();
        book(_h_pdg421_Dzero_pT_y->bin(i),     2, 1, y);
        book(_h_pdg411_Dplus_pT_y->bin(i),     3, 1, y);
        book(_h_pdg413_Dstarplus_pT_y->bin(i), 4, 1, y);
        book(_h_pdg431_Dsplus_pT_y->bin(i),    5, 1, y);
      }
      book(_h_pdg4122_Lambdac_pT ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Use PrimaryHadrons to avoid double counting and automatically remove the contributions from unstable?
      const UnstableParticles &ufs = apply<UnstableParticles> (event, "UFS");
      for (const Particle& p : ufs.particles() ) {

        // We're only interested in charm hadrons
        if (!p.isHadron() || !p.hasCharm()) continue;

        // Kinematic acceptance
        const double y = p.absrap(); ///< Double analysis efficiency with a "two-sided LHCb"
        const double pT = p.pT();

        // Fiducial acceptance of the measurements
        if (pT > 8.0*GeV || y < 2.0 || y > 4.5) continue;

        /// Experimental selection removes non-prompt charm hadrons: we ignore those from b decays
        if (p.fromBottom()) continue;

        switch (p.abspid()) {
        case 411:
          _h_pdg411_Dplus_pT_y->fill(y, pT/GeV);
          break;
        case 421:
          _h_pdg421_Dzero_pT_y->fill(y, pT/GeV);
          break;
        case 431:
          _h_pdg431_Dsplus_pT_y->fill(y, pT/GeV);
          break;
        case 413:
          _h_pdg413_Dstarplus_pT_y->fill(y, pT/GeV);
          break;
        case 4122:
          _h_pdg4122_Lambdac_pT->fill(pT/GeV);
          break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scale_factor = 0.5 * crossSection()/microbarn / sumOfWeights();
      scale(_h_pdg411_Dplus_pT_y, scale_factor);
      scale(_h_pdg421_Dzero_pT_y, scale_factor);
      scale(_h_pdg431_Dsplus_pT_y, scale_factor);
      scale(_h_pdg413_Dstarplus_pT_y, scale_factor);
      _h_pdg4122_Lambdac_pT->scaleW(scale_factor);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_pdg411_Dplus_pT_y;
    Histo1DGroupPtr _h_pdg421_Dzero_pT_y;
    Histo1DGroupPtr _h_pdg431_Dsplus_pT_y;
    Histo1DGroupPtr _h_pdg413_Dstarplus_pT_y;
    Histo1DPtr _h_pdg4122_Lambdac_pT;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCB_2013_I1218996);

}
