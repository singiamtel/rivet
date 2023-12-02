// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"

namespace Rivet {


  /// @brief CDF Run I Z \f$ p_\perp \f$ in Drell-Yan events
  ///
  /// @author Hendrik Hoeth
  class CDF_2000_S4155203 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2000_S4155203);


    /// @name Analysis methods
    /// @{

    void init() {
      // Set up projections
      DileptonFinder zfinder(91.2*GeV, 0.0, Cuts::abspid == PID::ELECTRON, Cuts::massIn(66*GeV, 116*GeV));
      declare(zfinder, "DileptonFinder");

      // Book histogram
      book(_hist_zpt ,1, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(e, "DileptonFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
        vetoEvent;
      }

      FourMomentum pZ = zfinder.bosons()[0].momentum();
      if (pZ.mass2() < 0) {
        MSG_DEBUG("Negative Z mass**2 = " << pZ.mass2()/GeV2 << "!");
        vetoEvent;
      }

      MSG_DEBUG("Dilepton mass = " << pZ.mass()/GeV << " GeV");
      MSG_DEBUG("Dilepton pT   = " << pZ.pT()/GeV << " GeV");
      _hist_zpt->fill(pZ.pT()/GeV);
    }


    void finalize() {
      scale(_hist_zpt, crossSection()/picobarn/sumOfWeights());
    }

    /// @}


  private:

    Histo1DPtr _hist_zpt;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2000_S4155203, CDF_2000_I505738);

}
