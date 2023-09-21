// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1989_I277658 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1989_I277658);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      offset = 0;
      if (isCompatibleWithSqrtS(14.0*GeV)) {
        offset = 1;
      }
      else if (isCompatibleWithSqrtS(22.0*GeV)) {
        offset = 2;
      }
      else if (isCompatibleWithSqrtS(34.8*GeV)) {
        offset = 3;
      }
      else if (isCompatibleWithSqrtS(43.6*GeV)) {
        offset = 4;
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      book(_histCh, 5, 1, offset);
      book(_histTotal, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (Ecm == "")  Ecm = _histTotal->bin(offset).xEdge();
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histCh->fill(cfs.size());
      _histTotal->fill(Ecm, cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histCh, 1.0/sumOfWeights());
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    BinnedHistoPtr<int> _histCh;
    BinnedProfilePtr<string> _histTotal;
    string Ecm = "";
    int offset;
    /// @}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1989_I277658);


}
