// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ATLAS inclusive jet pT spectrum, di-jet Mass and di-jet chi
  class ATLAS_2010_I871366: public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2010_I871366);

    /// Choice of jet algorithm
    enum Alg { AKT4=0, AKT6=1 };


    void init() {
      FinalState fs;
      declare(fs, "FinalState");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.6), "AntiKT06");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "AntiKT04");

      const vector<double> ybins{ 0.0, 0.3, 0.8, 1.2, 2.1, 2.8 };
      const vector<double> massBinsForChi{ 340, 520, 800, 1200 };


      size_t ptDsOffset(0), massDsOffset(10), chiDsOffset(20);
      for (size_t alg = 0; alg < 2; ++alg) {
        book(_pThistos[alg], ybins);
        for (auto& b : _pThistos[alg]->bins()) {
          book(b, b.index() + ptDsOffset, 1, 1);
        }
        ptDsOffset += 5;

        book(_massVsY[alg], ybins);
        for (auto& b : _massVsY[alg]->bins()) {
          book(b, b.index() + massDsOffset, 1, 1);
        }
        massDsOffset += 5;

        book(_chiVsMass[alg], massBinsForChi);
        for (auto& b : _chiVsMass[alg]->bins()) {
          book(b, b.index() + chiDsOffset, 1, 1);
        }
        chiDsOffset += 3;
      }
    }


    void analyze(const Event& evt) {
      Jets jetAr[2];
      jetAr[AKT6] = apply<FastJets>(evt, "AntiKT06").jetsByPt(Cuts::pT > 30*GeV);
      jetAr[AKT4] = apply<FastJets>(evt, "AntiKT04").jetsByPt(Cuts::pT > 30*GeV);

      // Identify the dijets
      for (size_t alg = 0; alg < 2; ++alg) {
        vector<FourMomentum> leadjets;
        for (const Jet& jet : jetAr[alg]) {
          const double pT = jet.pT();
          const double absy = jet.absrap();
          _pThistos[alg]->fill(absy, pT/GeV);

          if (absy < 2.8 && leadjets.size() < 2) {
            if (leadjets.empty() && pT < 60*GeV) continue;
            leadjets.push_back(jet.momentum());
          }
        }

        // Veto if no acceptable dijet found
        if (leadjets.size() < 2) {
          MSG_DEBUG("Could not find two suitable leading jets");
          continue;
        }

        const double rap1 = leadjets[0].rapidity();
        const double rap2 = leadjets[1].rapidity();
        const double mass = (leadjets[0] + leadjets[1]).mass();
        const double ymax = max(fabs(rap1), fabs(rap2));
        const double chi = exp(fabs(rap1 - rap2));
        if (fabs(rap1 + rap2) < 2.2) {
          _chiVsMass[alg]->fill(mass/GeV, chi);
        }
        _massVsY[alg]->fill(ymax, mass/GeV);

      }
    }


    void finalize() {
      const double sf = crossSectionPerEvent()/picobarn;
      // factor 0.5 needed because it is differential in dy and not d|y|
      scale(_pThistos, 0.5*sf);
      scale(_massVsY, sf);
      scale(_chiVsMass, sf);

      divByGroupWidth(_pThistos);
      divByGroupWidth(_massVsY);
      divByGroupWidth(_chiVsMass);
    }


  private:

    /// The inclusive pT spectrum for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _pThistos[2];

    /// The di-jet mass spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _massVsY[2];

    /// The di-jet chi distribution binned in mass for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _chiVsMass[2];

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2010_I871366, ATLAS_2010_S8817804);

}
