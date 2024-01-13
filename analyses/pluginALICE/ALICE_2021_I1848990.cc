// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief Measurement of beauty and charm production iat 5 TeV
  class ALICE_2021_I1848990 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1848990);

    void init() {

      const UnstableParticles fs(Cuts::absrap < 0.5);
      declare(fs, "fs");

      book(_h["npmD0"],    1, 1, 1); // non-prompt D0 production cross section
      book(_h["npmDplus"], 2, 1, 1); // non-prompt D+ production cross section
      book(_h["npmDS"],    3, 1, 1); // non-prompt DS+ production cross section
      book(_h["pmD0"],     4, 1, 1); // prompt D0 production cross section
      book(_h["pmDplus"],  5, 1, 1); // prompt D+ production cross section
      book(_h["pmDS"],     6, 1, 1); // prompt DS+ production cross section

      book(_e["npmD0_pmD0"],           7, 1, 1); // ratio between non-prompt D0 and prompt D0
      book(_e["npmDplus_pmDplus"],     8, 1, 1); // ratio between non-prompt D+ and prompt D+
      book(_e["npmDS_pmDS"],           9, 1, 1); // ratio between non-prompt DS+ and prompt DS+
      book(_e["pmDplus_pmD0"],        10, 1, 1); // ratio between prompt D+ and prompt D0
      book(_e["npmDplus_npmD0"],      11, 1, 1); // ratio between non-prompt D+ and non-prompt D0
      book(_e["pmDS_pmDplusD0sum"],   12, 1, 1); // ratio between prompt DS+ and sum of prompt D0 and prompt D+
      book(_e["npmDS_npmDplusD0sum"], 13, 1, 1); // ratio between non-prompt DS+ and sum of non-prompt D0 and non-prompt D+

      book(_h["pmDplusD0sum"],  "TMP/_pmDplusD0sum",  refData(12,1,1));
      book(_h["npmDplusD0sum"], "TMP/_npmDplusD0sum", refData(13,1,1));

      book(_e_frag["pmDS_pmDpD0"], 14, 1, 1); // c-quark fragmentation-fraction ratio
      book(_h_frag["pmDS"], "TMP/_fpmDS", refData<YODA::BinnedEstimate<string>>(14,1,1));
      book(_h_frag["pmDpD0"], "TMP/_fpmDpD0", refData<YODA::BinnedEstimate<string>>(14,1,1));

      book(_e_frag["npmDS_npmDpD0"], 15, 1, 1); // b-quark fragmentation-fraction ratio
      book(_h_frag["npmDS"], "TMP/_npmDS", refData<YODA::BinnedEstimate<string>>(15,1,1));
      book(_h_frag["npmDpD0"], "TMP/_npmDpD0", refData<YODA::BinnedEstimate<string>>(15,1,1));

    }

    void analyze(const Event& event) {

      const UnstableParticles& fs = apply<UnstableParticles>(event, "fs");

      for (const Particle& p : fs.particles()) {
        if (p.fromBottom()) {
          if (p.abspid() == PID::D0) {
            _h["npmD0"]->fill(p.pT()/GeV);
            _h["npmDplusD0sum"]->fill(p.pT()/GeV);
            if (p.pT() >= 2*GeV) _h_frag["npmDpD0"]->fill("PP --> D (Q=NON-PROMPT) + X"s);
          }
          else if (p.abspid() == PID::DPLUS) {
            _h["npmDplus"]->fill(p.pT()/GeV);
            _h["npmDplusD0sum"]->fill(p.pT()/GeV);
            if (p.pT() >= 2*GeV) _h_frag["npmDpD0"]->fill("PP --> D (Q=NON-PROMPT) + X"s);
          }
          else if (p.abspid() == PID::DSPLUS) {
            _h["npmDS"]->fill(p.pT()/GeV);
            if (p.pT() >= 2*GeV) _h_frag["npmDS"]->fill("PP --> D (Q=NON-PROMPT) + X"s);
          }
        }
        else {
          if (p.abspid() == PID::D0) {
            _h["pmD0"]->fill(p.pT()/GeV);
            _h["pmDplusD0sum"]->fill(p.pT()/GeV);
            if (p.pT() >= 1*GeV) _h_frag["pmDpD0"]->fill("PP --> D (Q=PROMPT) + X"s);
          }
          else if (p.abspid() == PID::DPLUS) {
            _h["pmDplus"]->fill(p.pT()/GeV);
            _h["pmDplusD0sum"]->fill(p.pT()/GeV);
            if (p.pT() >= 1*GeV) _h_frag["pmDpD0"]->fill("PP --> D (Q=PROMPT) + X"s);
          }
          else if (p.abspid() == PID::DSPLUS) {
            _h["pmDS"]->fill(p.pT()/GeV);
            if (p.pT() >= 1*GeV) _h_frag["pmDS"]->fill("PP --> D (Q=PROMPT) + X"s);
          }
        }
      }
    }

    void finalize() {

      const double sf = crossSection() / (microbarn*2*sumOfWeights());
      scale(_h, sf);
      scale(_h_frag, 0.5*sf);

      for (auto& item : _e) {
        size_t delim = item.first.find("_");
        YODA::Histo1D num = *_h[item.first.substr(0, delim)];
        if (!item.second->binning().isCompatible(num.binning())) {
          num = num.clone();
          num.rebinXTo(item.second->xEdges());
        }
        YODA::Histo1D den = *_h[item.first.substr(delim+1)];
        if (!item.second->binning().isCompatible(den.binning())) {
          den = den.clone();
          den.rebinXTo(item.second->xEdges());
        }
        divide(num, den, item.second);
      }

      for (auto& item : _e_frag) {
        size_t delim = item.first.find("_");
        divide(_h_frag[item.first.substr(0, delim)],
               _h_frag[item.first.substr(delim+1)],
               item.second);
      }
    }

  private:

    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<string>> _h_frag;

    map<string, Estimate1DPtr> _e;
    map<string, BinnedEstimatePtr<string>> _e_frag;

  };

  RIVET_DECLARE_PLUGIN(ALICE_2021_I1848990);

}
