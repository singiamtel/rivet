// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2013_I1234228 : public Analysis {
  public:


    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1234228);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      DileptonFinder zfinder(91.2*GeV, 0.1, cuts && Cuts::abspid == PID::ELECTRON, Cuts::massIn(116*GeV, 1500*GeV));
      declare(zfinder, "DileptonFinder");

      book(_hist_mll, 1, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(event, "DileptonFinder");

      if (zfinder.bosons().size() != 1)  vetoEvent;

      double mass = zfinder.bosons()[0].mass();
      _hist_mll->fill(mass);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/sumOfWeights();
      scale(_hist_mll, sf);
    }

    /// @}

  private:


    /// @name Histograms
    /// @{
    Histo1DPtr _hist_mll;
    /// @}

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1234228);

}
