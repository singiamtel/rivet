// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi, and psi(2s) production at 5.02 and 13 TeV
  class ALICE_2017_I1511865 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2017_I1511865);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::pid==443 || Cuts::pid==100443), "UFS");

      for (double eVal : allowedEnergies()) {
        const string en = toString(int(eVal/MeV));
        if (isCompatibleWithSqrtS(eVal))  _sqs = en;

        size_t offset = (en == "5020"s)? 6 : 0;
        book(_h[en+"JPsi_pT"], 1+offset,1,1);
        book(_h[en+"JPsi_y"],  2+offset,1,1);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
      book(_h["Psi2S_pT"], 3,1,1);
      book(_h["Psi2S_y"],  4,1,1);
      book(_h["JPsi_pT2"], "TMP/JPsi_pY", refData(5,1,1));
      book(_h["JPsi_y2"],  "TMP/JPsi_y",  refData(6,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        // rapidity cut
        const double absrap = p.absrap();
        if (absrap<2.5 || absrap>4) continue;
        const double xp = p.perp();
        // J/Psi
        if (p.pid()==443) {
          if (_sqs != "5020"s) {
            if (xp>30.) continue;
            _h[_sqs+"JPsi_pT"]->fill(xp);
            _h[_sqs+"JPsi_y"]->fill(absrap);
            if (xp<=16.) {
              _h["JPsi_pT2"]->fill(xp);
              _h["JPsi_y2"]->fill(absrap);
            }
          }
          else {
            if (xp>12.) continue;
            _h[_sqs+"JPsi_pT"]->fill(xp);
            _h[_sqs+"JPsi_y"]->fill(absrap);
          }
        }
        // psi(2S)
        else if (_sqs != "5020"s) {
          if (xp>16.) continue;
          _h["Psi2S_pT"]->fill(xp);
          _h["Psi2S_y"]->fill(absrap);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor 1/2 due folding +/- rap
      const double fact = 0.5*crossSection()/nanobarn/sumOfWeights();
      // factor 1.5 for rapidity range 2.5-4
      for (auto& item : _h) {
        if (item.first.find("_pT") != string::npos) {
          scale(item.second, fact/1.5);
        }
        else {
          scale(item.second, fact);
        }
      }

      Estimate1DPtr tmp;
      book(tmp,5,1,1);
      divide(_h["Psi2S_pT"], _h["JPsi_pT2"], tmp);
      book(tmp,6,1,1);
      divide(_h["Psi2S_y"], _h["JPsi_y2"], tmp);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string,Histo1DPtr> _h;
    Histo1DPtr _h_JPsi_pT[2], _h_JPsi_y[2];
    Histo1DPtr _h_Psi2S_pT,_h_Psi2S_y;
    Histo1DPtr _h_JPsi_pT2,_h_JPsi_y2;

    string _sqs = "";
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2017_I1511865);

}
