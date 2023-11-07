// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class ATLAS_2012_I1082936 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2012_I1082936);

    /// @}


  public:

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");

      FastJets fj04(fs,  JetAlg::ANTIKT, 0.4);
      fj04.useInvisibles();
      declare(fj04, "AntiKT04");

      FastJets fj06(fs,  JetAlg::ANTIKT, 0.6);
      fj06.useInvisibles();
      declare(fj06, "AntiKT06");


      // Histogram booking copied from the previous analysis
      const vector<double> ybins{ 0.0, 0.3, 0.8, 1.2, 2.1, 2.8, 3.6, 4.4 };
      const vector<double> ystarbins{ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.4};

      size_t ptDsOffset(0), massDsOffset(2);
      for (size_t alg = 0; alg < 2; ++alg) {
        book(_pThistos[alg], ybins);
        for (auto& b : _pThistos[alg]->bins()) {
          book(b, 1+ptDsOffset, 1, b.index());
        }
        ++ptDsOffset;

        book(_mass[alg], ystarbins);
        for (auto& b : _mass[alg]->bins()) {
          book(b, 1+massDsOffset, 1, b.index());
        }
        ++massDsOffset;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Jets jetAr[2];
      jetAr[AKT6] = apply<FastJets>(event, "AntiKT06").jetsByPt(20*GeV);
      jetAr[AKT4] = apply<FastJets>(event, "AntiKT04").jetsByPt(20*GeV);

      // Loop over jet "radii" used in analysis
      for (size_t alg = 0; alg < 2; ++alg) {
        // Identify dijets
        vector<FourMomentum> leadjets;
        for (const Jet& jet : jetAr[alg]) {
          const double pT = jet.pT();
          const double absy = jet.absrap();
          _pThistos[alg]->fill(absy, pT/GeV);

          if (absy < 4.4 && leadjets.size() < 2) {
            if (leadjets.empty() && pT < 30*GeV) continue;
            leadjets.push_back(jet.momentum());
          }
        }
        // Make sure we have a leading jet with pT >30 GeV and a second to leading jet with pT>20 GeV
        if (leadjets.size() < 2) {
          MSG_DEBUG("Could not find two suitable leading jets");
          continue;
        }

        const double y1 = leadjets[0].rapidity();
        const double y2 = leadjets[1].rapidity();
        const double ystar = fabs(y1-y2)/2.;
        const double m = (leadjets[0] + leadjets[1]).mass();
        // Fill mass histogram
        _mass[alg]->fill(ystar, m/TeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor 0.5 needed because it is differential in dy and not d|y|
      scale(_pThistos, 0.5*crossSectionPerEvent()/picobarn);
      scale(_mass, crossSectionPerEvent()/picobarn);
      divByGroupWidth(_pThistos);
      divByGroupWidth(_mass);
    }

    /// @}


  private:

    // Data members like post-cuts event weight counters go here

    enum Alg { AKT4=0, AKT6=1 };

  private:

    /// The inclusive pT spectrum for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _pThistos[2];

    /// The di-jet mass spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _mass[2];
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2012_I1082936);

}
