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

      size_t ih = 0;
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal/MeV));
        if (isCompatibleWithSqrtS(eVal)) {
          _sqs = en;
          _ie = ih+1;
        }
        book(_h[en], 5, 1, ++ih);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
      book(_p, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (Ecm == "")  Ecm = _p->bin(_ie).xEdge();
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _h[_sqs]->fill(cfs.size());
      _p->fill(Ecm, cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1.0/sumOfWeights());
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    map<string,BinnedHistoPtr<int>> _h;
    BinnedProfilePtr<string> _p;
    string Ecm = "", _sqs = "";
    size_t _ie = 0;
    /// @}
  };


  RIVET_DECLARE_PLUGIN(TASSO_1989_I277658);


}
