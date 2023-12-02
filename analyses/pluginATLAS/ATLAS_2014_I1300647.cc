// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, DileptonFinder...

namespace Rivet {


  class ATLAS_2014_I1300647 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1300647);


  public:

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      DileptonFinder zfinder_dressed_el(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_dressed_el, "DileptonFinder_dressed_el");

      DileptonFinder zfinder_bare_el(91.2*GeV, 0.0, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV && Cuts::abspid == PID::ELECTRON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_bare_el,	"DileptonFinder_bare_el");

      DileptonFinder zfinder_dressed_mu(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_dressed_mu, "DileptonFinder_dressed_mu");

      DileptonFinder zfinder_bare_mu(91.2*GeV, 0.0, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV && Cuts::abspid == PID::MUON, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder_bare_mu,	"DileptonFinder_bare_mu");

      // Book histograms
      book(_hist_zpt_el_dressed ,1, 1, 1);  // electron "dressed"
      book(_hist_zpt_mu_dressed ,1, 1, 2);  // muon "dressed"
      book(_hist_zpt_el_bare    ,1, 2, 1);  // electron "bare"
      book(_hist_zpt_mu_bare    ,1, 2, 2);  // muon "bare"

      //double-differential plots
      book(_h_zpt_el_mu_dressed, {0., 1., 2., 2.4}, {"d03-x01-y02", "d03-x01-y04", "d03-x01-y06"});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DileptonFinder& zfinder_dressed_el = apply<DileptonFinder>(event, "DileptonFinder_dressed_el");
      const DileptonFinder& zfinder_bare_el    = apply<DileptonFinder>(event, "DileptonFinder_bare_el");
      const DileptonFinder& zfinder_dressed_mu = apply<DileptonFinder>(event, "DileptonFinder_dressed_mu");
      const DileptonFinder& zfinder_bare_mu    = apply<DileptonFinder>(event, "DileptonFinder_bare_mu");

      FillPlots1d(zfinder_dressed_el, _hist_zpt_el_dressed);

      FillPlots1d(zfinder_bare_el,    _hist_zpt_el_bare);

      FillPlots1d(zfinder_dressed_mu, _hist_zpt_mu_dressed);

      FillPlots1d(zfinder_bare_mu,    _hist_zpt_mu_bare);

      FillPlots3d(zfinder_dressed_el, _h_zpt_el_mu_dressed);
      FillPlots3d(zfinder_dressed_mu, _h_zpt_el_mu_dressed);

    }

    void FillPlots1d(const DileptonFinder& zfinder, Histo1DPtr hist) {
      if(zfinder.bosons().size() != 1) return;
      const FourMomentum pZ = zfinder.bosons()[0].momentum();
      hist->fill(pZ.pT()/GeV);
      return;
    }

    void FillPlots3d(const DileptonFinder& zfinder, Histo1DGroupPtr& binnedHist) {
      if(zfinder.bosons().size() != 1) return;
      const FourMomentum pZ = zfinder.bosons()[0].momentum();
      binnedHist->fill(pZ.rapidity(), pZ.pT()/GeV);
      return;
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_hist_zpt_el_dressed);
      normalize(_hist_zpt_el_bare);

      normalize(_hist_zpt_mu_dressed);
      normalize(_hist_zpt_mu_bare);

      normalize(_h_zpt_el_mu_dressed);

    }

    /// @}


  private:

    // Data members like post-cuts event weight counters go here


  private:

    /// @name Histograms
    /// @{
    Histo1DGroupPtr _h_zpt_el_mu_dressed;


    Histo1DPtr _hist_zpt_el_dressed;
    Histo1DPtr _hist_zpt_el_bare;
    Histo1DPtr _hist_zpt_mu_dressed;
    Histo1DPtr _hist_zpt_mu_bare;

    /// @}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1300647);
}
