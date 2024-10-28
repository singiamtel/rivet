// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// UE charged particles vs. leading jet
  class CMS_2011_I916908 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I916908);


    void init() {
      const ChargedFinalState cfs(Cuts::abseta < 2.0 && Cuts::pT >= 500*MeV);
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet(Cuts::abseta < 2.5 && Cuts::pT >= 500*MeV);
      const FastJets jetpro(cfsforjet, JetAlg::SISCONE, 0.5);
      declare(jetpro, "Jets");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        if (en == "7000"s) {
          book(_p[en+"Nch_vs_pT"], 1, 1, 1); // Nch vs. pT_max
          book(_p[en+"Sum_vs_pT"], 2, 1, 1); // sum(pT) vs. pT_max
          book(_h[en+"pT3_Nch"],   5, 1, 1);   // transverse Nch,     pT_max > 3GeV
          book(_h[en+"pT3_Sum"],   6, 1, 1);   // transverse sum(pT), pT_max > 3GeV
          book(_h[en+"pT3_pT"],    7, 1, 1);   // transverse pT,      pT_max > 3GeV
          book(_h[en+"pT20_Nch"],  8, 1, 1);   // transverse Nch,     pT_max > 20GeV
          book(_h[en+"pT20_Sum"],  9, 1, 1);   // transverse sum(pT), pT_max > 20GeV
          book(_h[en+"pT20_pT"],  10, 1, 1);  // transverse pT,      pT_max > 20GeV
        }
        else if (en == "900"s) {
          book(_p[en+"Nch_vs_pT"], 3, 1, 1); // Nch vs. pT_max
          book(_p[en+"Sum_vs_pT"], 4, 1, 1); // sum(pT) vs. pT_max
          book(_h[en+"pT3_Nch"],  11, 1, 1);  // transverse Nch,     pT_max > 3GeV
          book(_h[en+"pT3_Sum"],  12, 1, 1);  // transverse sum(pT), pT_max > 3GeV
          book(_h[en+"pT3_pT"],   13, 1, 1);  // transverse pT,      pT_max > 3GeV
        }
        book(_c[en+"pT3"],  "TMP/nch_tot_pT3"+en);
        book(_c[en+"pT20"], "TMP/nch_tot_pT20"+en);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_lead;
      for (const Jet& j : apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 1.0*GeV && Cuts::abseta < 2.0)) {
        p_lead = j.momentum();
        break;
      }
      if (p_lead.isZero()) vetoEvent;
      const double philead = p_lead.phi();
      const double pTlead  = p_lead.pT();

      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();

      int nTransverse = 0;
      double ptSumTransverse = 0.;
      for (const Particle& p : particles) {
        double dphi = deltaPhi(philead, p.phi());
        if (dphi>PI/3. && dphi<PI*2./3.) {   // Transverse region
          ++nTransverse;

          const double pT = p.pT()/GeV;
          ptSumTransverse += pT;

          if (pTlead > 3*GeV) _h[_sqs+"pT3_pT"]->fill(pT/GeV);
          if (_sqs == "7000"s && pTlead > 20*GeV) _h[_sqs+"pT20_pT"]->fill(pT/GeV);
        }
      }

      const double area = 8./3. * PI;
      _p[_sqs+"Nch_vs_pT"]->fill(pTlead/GeV, 1./area*nTransverse);
      _p[_sqs+"Sum_vs_pT"]->fill(pTlead/GeV, 1./area*ptSumTransverse);
      if (pTlead > 3.0*GeV) {
        _h[_sqs+"pT3_Nch"]->fill(nTransverse);
        _h[_sqs+"pT3_Sum"]->fill(ptSumTransverse/GeV);
        _c[_sqs+"pT3"]->fill(nTransverse);
      }
      if (_sqs == "7000"s && pTlead > 20*GeV) {
        _h[_sqs+"pT20_Nch"]->fill(nTransverse);
        _h[_sqs+"pT20_Sum"]->fill(ptSumTransverse/GeV);
        _c[_sqs+"pT20"]->fill(nTransverse);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal));
        if (_h[en+"pT3_pT"]->sumW()) {
          normalize(_h[en+"pT3_pT"], dbl(*_c[en+"pT3"]) / _h[en+"pT3_Nch"]->sumW());
        }
        if (en == "7000"s && _h[en+"pT20_pT"]->sumW()) {
          normalize(_h[en+"pT20_pT"], dbl(*_c[en+"pT20"]) / _h[en+"pT20_Nch"]->sumW());
        }
      }
      for (auto& item : _h) {
        if (item.first.find("_pT") != string::npos)  continue;
        normalize(item.second);
      }

    }



  private:

    /// @{
    map<string,CounterPtr> _c;
    map<string,Profile1DPtr> _p;
    map<string,Histo1DPtr> _h;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_I916908, CMS_2011_S9120041);

}
