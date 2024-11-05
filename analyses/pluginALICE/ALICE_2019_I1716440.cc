// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief D meson production at 5.02 and 7 TeV
  class ALICE_2019_I1716440 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2019_I1716440);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==411 || Cuts::abspid==421 || Cuts::abspid==431 || Cuts::abspid==413), "UFS");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal/GeV));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;
        bool is7TeV(en == "7000"s);

        size_t offset = is7TeV? 10 : 6;
        for (size_t ix=0; ix<4; ++ix) {
          const string mode = en+toString(ix+1);
          book(_h[mode+"ratio_num"], "TMP/h_ratio_num_"+mode, refData(offset+ix,1,1));
          book(_h[mode+"ratio_den"], "TMP/h_ratio_den_"+mode, refData(offset+ix,1,1));
          book(_h[mode+"energy"],    "h_energy_"+mode,        refData(offset+ix,1,1));
          if (!is7TeV)  book(_h[mode+"prompt"], 1+ix, 1, 1);
        }
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
      book(_h["incl"], 5, 1, 1);
      book(_total, 18, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        // no mixing and |y|<0.5
        if (p.children().size()==1 || p.absrap()>0.5) continue;
        unsigned int imeson=0;
        if (p.abspid()==411)       imeson=1;
        else if (p.abspid()==413)  imeson=2;
        else if (p.abspid()==431)  imeson=3;
        const double pT = p.perp();
        if (p.fromBottom()) {
          if (_sqs == "5020"s && imeson==0) _h["incl"]->fill(pT);
          continue;
        }
        // prompt at 5.02 TeV
        if (_sqs == "7000"s) {
          _h[_sqs+toString(imeson)+"prompt"]->fill(pT);
          _total->fill(_edges[imeson]);
        }
        if (imeson==0) {
          _h[_sqs+"0ratio_den"]->fill(pT/GeV);
          _h[_sqs+"1ratio_den"]->fill(pT/GeV);
          _h[_sqs+"2ratio_den"]->fill(pT/GeV);
          if (_sqs == "5020"s) _h["incl"]->fill(pT/GeV);
        }
        else if (imeson==1) {
          _h[_sqs+"0ratio_num"]->fill(pT/GeV);
          _h[_sqs+"3ratio_den"]->fill(pT/GeV);
        }
        else if (imeson==2) {
          _h[_sqs+"1ratio_num"]->fill(pT/GeV);
        }
        else if (imeson==3) {
          _h[_sqs+"2ratio_num"]->fill(pT/GeV);
          _h[_sqs+"3ratio_num"]->fill(pT/GeV);
        }
        _h[_sqs+toString(imeson)+"energy"]->fill(pT/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double factor = crossSection()/microbarn/sumOfWeights();
      scale(_total, factor);
      scale(_h,factor);

      Estimate1DPtr tmp;
      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal/GeV));
        for (size_t ix=0; ix<4; ++ix) {
          const string mode = en+toString(ix+1);
          bool is7TeV(en == "7000"s);
          book(tmp, 6+4*is7TeV+ix, 1, 1);
          divide(_h[mode+"ratio_num"], _h[mode+"ratio_den"], tmp);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h;
    BinnedHistoPtr<string> _total;
    vector<string> _edges = { "P P --> D0 (Q=PROMPT) X",
                              "P P --> D+ (Q=PROMPT) X",
                              "P P --> D*+ (Q=PROMPT) X",
                              "P P --> D/s+ (Q=PROMPT) X" };
    string _sqs = "";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2019_I1716440);

}
