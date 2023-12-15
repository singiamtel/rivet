// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Particle.fhh"

namespace Rivet {


  /// Low-mass Drell-Yan differential cross section at 7 TeV
  class ATLAS_2014_I1288706 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1288706);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Set up projections
      DileptonFinder zfinder_ext_dressed_mu(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 6.0*GeV && Cuts::abspid == PID::MUON, Cuts::massIn(12.0*GeV, 66.0*GeV));
      declare(zfinder_ext_dressed_mu, "DileptonFinder_ext_dressed_mu");

      DileptonFinder zfinder_dressed_mu(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 12*GeV && Cuts::abspid == PID::MUON, Cuts::massIn(26.0*GeV, 66.0*GeV));
      declare(zfinder_dressed_mu, "DileptonFinder_dressed_mu");

      DileptonFinder zfinder_dressed_el(91.2*GeV, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 12*GeV && Cuts::abspid == PID::ELECTRON, Cuts::massIn(26.0*GeV, 66.0*GeV));
      declare(zfinder_dressed_el, "DileptonFinder_dressed_el");

      book(_hist_ext_mu_dressed, 1, 1, 1); // muon, dressed level, extended phase space
      book(_hist_mu_dressed, 2, 1, 1); // muon, dressed level, nominal phase space
      book(_hist_el_dressed, 2, 1, 2); // electron, dressed level, nominal phase space
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const DileptonFinder& zfinder_ext_dressed_mu = apply<DileptonFinder>(event, "DileptonFinder_ext_dressed_mu");
      const DileptonFinder& zfinder_dressed_mu     = apply<DileptonFinder>(event, "DileptonFinder_dressed_mu"    );
      const DileptonFinder& zfinder_dressed_el     = apply<DileptonFinder>(event, "DileptonFinder_dressed_el"    );

      fillPlots(zfinder_ext_dressed_mu, _hist_ext_mu_dressed, 9*GeV);
      fillPlots(zfinder_dressed_mu,     _hist_mu_dressed,    15*GeV);
      fillPlots(zfinder_dressed_el,     _hist_el_dressed,    15*GeV);
    }


    /// Helper function to fill the histogram if a Z is found with the required lepton cuts
    void fillPlots(const DileptonFinder& zfinder, Histo1DPtr hist, double leading_pT) {
      if (zfinder.bosons().size() != 1) return;
      const FourMomentum el1 = zfinder.leptons()[0];
      const FourMomentum el2 = zfinder.leptons()[1];
      if (el1.pT() < leading_pT && el2.pT() < leading_pT) return;
      const double mass = zfinder.bosons()[0].mass();
      hist->fill(mass/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_ext_mu_dressed, crossSection()/sumOfWeights());
      scale(_hist_mu_dressed,     crossSection()/sumOfWeights());
      scale(_hist_el_dressed,     crossSection()/sumOfWeights());
    }

    /// @}


  private:

    /// Histograms
    Histo1DPtr _hist_ext_mu_dressed, _hist_mu_dressed, _hist_el_dressed;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1288706);

}
