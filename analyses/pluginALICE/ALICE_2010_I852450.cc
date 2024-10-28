// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2010_I852450 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2010_I852450);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs05(Cuts::abseta < 0.5);
      ChargedFinalState cfs10(Cuts::abseta < 1.0);
      ChargedFinalState cfs13(Cuts::abseta < 1.3);
      declare(cfs05, "CFS05");
      declare(cfs10, "CFS10");
      declare(cfs13, "CFS13");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        size_t offset = (en == "2360"s)? 6 : 0;

        book(_h[en+"dN_dNch_05"], 11+offset, 1, 1);
        book(_h[en+"dN_dNch_10"], 12+offset, 1, 1);
        book(_h[en+"dN_dNch_13"], 13+offset, 1, 1);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& charged_05 = apply<ChargedFinalState>(event, "CFS05");
      const ChargedFinalState& charged_10 = apply<ChargedFinalState>(event, "CFS10");
      const ChargedFinalState& charged_13 = apply<ChargedFinalState>(event, "CFS13");

      _h[_sqs+"dN_dNch_05"]->fill(charged_05.size());
      _h[_sqs+"dN_dNch_10"]->fill(charged_10.size());
      _h[_sqs+"dN_dNch_13"]->fill(charged_13.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h;

    string _sqs = "";
    /// @}


  };


  RIVET_DECLARE_ALIASED_PLUGIN(ALICE_2010_I852450, ALICE_2010_S8624100);

}
