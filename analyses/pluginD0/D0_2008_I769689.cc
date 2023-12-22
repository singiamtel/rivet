// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief D0 Run II Z \f$ p_\perp \f$ differential cross-section shape
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2008_I769689 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2008_I769689);


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      DileptonFinder zfinder(91.2*GeV, 0.2, Cuts::abspid == PID::ELECTRON, Cuts::massIn(40*GeV, 200*GeV));
      declare(zfinder, "DileptonFinder");

      book(_h_ZpT         ,1, 1, 1);
      book(_h_forward_ZpT ,3, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(e, "DileptonFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("No unique lepton pair found.");
        vetoEvent;
      }
      const double yZ = fabs(zfinder.bosons()[0].rapidity());
      const double pTZ = zfinder.bosons()[0].pT();
      _h_ZpT->fill(pTZ);
      if (yZ > 2) _h_forward_ZpT->fill(pTZ);
    }


    // Finalize
    void finalize() {
      normalize(_h_ZpT);
      normalize(_h_forward_ZpT);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_ZpT, _h_forward_ZpT;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2008_I769689, D0_2008_S7554427);

}
