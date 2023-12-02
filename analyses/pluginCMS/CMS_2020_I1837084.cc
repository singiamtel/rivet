// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// Measurements of differential Z-boson -> ll and vv production cross-sections in 13 TeV pp collisions
  class CMS_2020_I1837084 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2020_I1837084);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      DileptonFinder zmmFind(91.2*GeV, 0.1, Cuts::abspid == PID::MUON, Cuts::massIn(76.1876*GeV, 106.1876*GeV));
      declare(zmmFind, "ZmmFind");

      // Book histograms
      book(_h_Z_pt,      12, 1, 1);
      book(_e_Z_pt_norm, 13, 1, 1);
      book(_h_Z_pt_norm, "_h_Z_pt_norm", refData(13, 1, 1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& zmms = apply<DileptonFinder>(event, "ZmmFind").bosons();

      if (zmms.size() == 1 && zmms[0].pT() > 200*GeV) {
        _h_Z_pt     ->fill(min(zmms[0].pT()/GeV, 1499.999));
        _h_Z_pt_norm->fill(min(zmms[0].pT()/GeV, 1499.999));
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/femtobarn/sumOfWeights() : 1.0;

      scale(_h_Z_pt, norm);

      normalize(_h_Z_pt_norm);
      barchart(_h_Z_pt_norm, _e_Z_pt_norm);

    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Z_pt, _h_Z_pt_norm;
    Estimate1DPtr _e_Z_pt_norm;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2020_I1837084);

}
