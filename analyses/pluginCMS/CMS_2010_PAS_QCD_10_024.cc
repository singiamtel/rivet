// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  /// @brief Pseudorapidity distributions of charged particles at sqrt{s} = 0.9 and 7 TeV
  class CMS_2010_PAS_QCD_10_024 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2010_PAS_QCD_10_024);


    void init() {

      declare(ChargedFinalState(Cuts::abseta < 0.8 && Cuts::pT >= 0.5*GeV), "CFS_08_05");
      declare(ChargedFinalState(Cuts::abseta < 0.8 && Cuts::pT >= 1.0*GeV), "CFS_08_10");
      declare(ChargedFinalState(Cuts::abseta < 2.4 && Cuts::pT >= 0.5*GeV), "CFS_24_05");
      declare(ChargedFinalState(Cuts::abseta < 2.4 && Cuts::pT >= 1.0*GeV), "CFS_24_10");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        size_t offset = (en == "900"s)? 4 : 0;

        book(_h[en+"dNch_deta_pt05_eta08"], 1+offset, 1, 1);
        book(_h[en+"dNch_deta_pt10_eta08"], 2+offset, 1, 1);
        book(_h[en+"dNch_deta_pt05_eta24"], 3+offset, 1, 1);
        book(_h[en+"dNch_deta_pt10_eta24"], 4+offset, 1, 1);

        book(_c[en+"pt05_eta08"], "_pt05_eta08_"+en);
        book(_c[en+"pt10_eta08"], "_pt10_eta08_"+en);
        book(_c[en+"pt05_eta24"], "_pt05_eta24_"+en);
        book(_c[en+"pt10_eta24"], "_pt10_eta24_"+en);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    void analyze(const Event& event) {

      const ChargedFinalState& cfs_08_05 = apply<ChargedFinalState>(event, "CFS_08_05");
      const ChargedFinalState& cfs_08_10 = apply<ChargedFinalState>(event, "CFS_08_10");
      const ChargedFinalState& cfs_24_05 = apply<ChargedFinalState>(event, "CFS_24_05");
      const ChargedFinalState& cfs_24_10 = apply<ChargedFinalState>(event, "CFS_24_10");

      // Plot distributions
      if (!cfs_08_05.particles().empty()) _c[_sqs+"pt05_eta08"]->fill();
      if (!cfs_24_05.particles().empty()) _c[_sqs+"pt05_eta24"]->fill();
      for (const Particle& p : cfs_24_05.particles()) {
        _h[_sqs+"dNch_deta_pt05_eta24"]->fill(p.eta());
        if (!cfs_08_05.particles().empty()) {
          _h[_sqs+"dNch_deta_pt05_eta08"]->fill(p.eta());
        }
      }
      if (!cfs_08_10.particles().empty()) _c[_sqs+"pt10_eta08"]->fill();
      if (!cfs_24_10.particles().empty()) _c[_sqs+"pt10_eta24"]->fill();
      for (const Particle& p : cfs_24_10.particles()) {
        _h[_sqs+"dNch_deta_pt10_eta24"]->fill(p.eta());
        if (!cfs_08_10.particles().empty()) {
          _h[_sqs+"dNch_deta_pt10_eta08"]->fill(p.eta());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (dbl(*_c[en+"pt05_eta08"])) scale(_h[en+"dNch_deta_pt05_eta08"], 1./ *_c[en+"pt05_eta08"]);
        if (dbl(*_c[en+"pt10_eta08"])) scale(_h[en+"dNch_deta_pt10_eta08"], 1./ *_c[en+"pt10_eta08"]);
        if (dbl(*_c[en+"pt05_eta24"])) scale(_h[en+"dNch_deta_pt05_eta24"], 1./ *_c[en+"pt05_eta24"]);
        if (dbl(*_c[en+"pt10_eta24"])) scale(_h[en+"dNch_deta_pt10_eta24"], 1./ *_c[en+"pt10_eta24"]);
      }
    }


  private:

    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    string _sqs = "";
  };


  RIVET_DECLARE_PLUGIN(CMS_2010_PAS_QCD_10_024);

}
