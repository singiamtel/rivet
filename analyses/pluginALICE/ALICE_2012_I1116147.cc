//-*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class ALICE_2012_I1116147 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2012_I1116147);


    /// Initialise projections and histograms
    void init() {

      const UnstableParticles ufs(Cuts::absrap < RAPMAX);
      declare(ufs, "UFS");

      // Book histos
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        bool is900GeV(en == "900"s);
        book(_h[en+"pi0"], 1 + is900GeV,1,1);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
      book(_h["eta"], 3,1,1);
      // Temporary plots with the binning of _h_etaToPion to construct the eta/pi0 ratio
      book(_h["temp_eta"],  "TMP/h_eta" , refData(3,1,1));
      book(_h["temp_pion"], "TMP/h_pion", refData(3,1,1));
    }


    /// Per-event analysis
    void analyze(const Event& event) {

      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double normfactor = TWOPI*p.pT()/GeV*2*RAPMAX;
        if (p.pid() == 111) {
          // Neutral pion; ALICE corrects for pi0 feed-down from K_0_s and Lambda
          if (p.hasAncestorWith(Cuts::pid == 310)  ||
              p.hasAncestorWith(Cuts::pid == 3122) ||
              p.hasAncestorWith(Cuts::pid == -3122)) continue; //< K_0_s, Lambda, Anti-Lambda
          _h[_sqs+"pi0"]->fill(p.pT()/GeV, 1.0/normfactor);
          if (_sqs == "7000"s)  _h["temp_pion"]->fill(p.pT()/GeV);
        }
        else if (p.pid() == 221 && _sqs == "7000"s) {
          // eta meson (only for 7 TeV)
          _h["eta"]->fill(p.pT()/GeV, 1.0/normfactor);
          _h["temp_eta"]->fill(p.pT()/GeV);
        }
      }
    }


    /// Normalize histos and construct ratio
    void finalize() {
      const double sf = crossSection()/microbarn/sumOfWeights();
      scale(_h, sf);

      // first divide hists with binned axis
      Estimate1DPtr tmp_ratio;
      book(tmp_ratio, "TMP/ratio", _h["eta"]->xEdges());
      divide(_h["temp_eta"], _h["temp_pion"], tmp_ratio);
      // now convert to strings edges from ref data
      BinnedEstimatePtr<string> ratio;
      book(ratio, 4, 1, 1);
      for (const auto& b : tmp_ratio->bins()) {
        ratio->bin(b.index()).set(b.val(), b.err());
      }
    }


  private:

    const double RAPMAX = 0.8;
    string _sqs = "";

    map<string,Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(ALICE_2012_I1116147);

}
