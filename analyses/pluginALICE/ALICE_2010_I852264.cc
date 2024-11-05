// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ALICE_2010_I852264 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2010_I852264);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs(Cuts::abseta < 1.0);
      declare(cfs, "CFS");

      size_t ih = 0;
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        book(_h[en+"dN_deta"], 4+ih, 1, 1);
        if (ih == 2)  book(_h["dN_dNch"], 3, 1, 1);
        book(_c[en], "Nevt_after_cuts_"+en);
        ++ih;
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      if (charged.size() < 1)  vetoEvent;
      _c[_sqs]->fill();
      if (_sqs == "7000"s)  _h["dN_dNch"]->fill(charged.size());

      for (const Particle& p : charged.particles()) {
        const double eta = p.eta();
        _h[_sqs+"dN_deta"]->fill(eta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (_c[en]->sumW() == 0)  continue;
        scale(_h[en+"dN_deta"], 1.0/ *_c[en]);
      }
      normalize(_h["dN_dNch"]);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h;
    map<string, CounterPtr> _c;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ALICE_2010_I852264, ALICE_2010_S8625980);

}
