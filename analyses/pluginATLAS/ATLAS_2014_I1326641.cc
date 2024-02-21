// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2014_I1326641 : public Analysis {
  public:

    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1326641);

    /// @}


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      FastJets fj04(fs, JetAlg::ANTIKT, 0.4);
      fj04.useInvisibles();
      declare(fj04, "AntiKT04");

      FastJets fj06(fs, JetAlg::ANTIKT, 0.6);
      fj06.useInvisibles();
      declare(fj06, "AntiKT06");

      const vector<double> ystarBins{ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };

      size_t massDsOffset(0);
      for (size_t alg = 0; alg < 2; ++alg) {
        book(h_trijet_Mass[alg], ystarBins);
        for (auto& b : h_trijet_Mass[alg]->bins()) {
          book(b, 1+massDsOffset, 1, 1);
          ++massDsOffset;
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Jets jetAr[2];
      jetAr[AKT4] = apply<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT > 50*GeV);
      jetAr[AKT6] = apply<FastJets>(event, "AntiKT06").jetsByPt(Cuts::pT > 50*GeV);

      const size_t nJets = 3;
      double ptCut[nJets] = { 150., 100., 50.};

      // Loop over jet "radii" used in analysis
      for (size_t alg = 0; alg < 2; ++alg) {
        // Identify 3jets
        vector<FourMomentum> leadJets;
        for (const Jet& jet : jetAr[alg]) {
          if (jet.absrap() < 3.0 && leadJets.size() < nJets) {
            int filledJets = leadJets.size();
            if (jet.pT() < ptCut[filledJets])  continue;
            leadJets.push_back(jet.momentum());
          }
        }

        if (leadJets.size() < nJets) {
          MSG_DEBUG("Could not find three suitable leading jets");
          continue;
        }

        const double y1 = leadJets[0].rapidity();
        const double y2 = leadJets[1].rapidity();
        const double y3 = leadJets[2].rapidity();

        const double yStar = fabs(y1-y2) + fabs(y2-y3) + fabs(y1-y3);
        const double m = (leadJets[0] + leadJets[1] + leadJets[2]).mass();
        h_trijet_Mass[alg]->fill(yStar, m);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf( 2.0 * crossSection()/picobarn / sumOfWeights());
      scale(h_trijet_Mass, sf);

    }

    /// @}


  private:

    // Data members like post-cuts event weight counters go here
    enum Alg { AKT4=0, AKT6=1 };

  private:

    // The 3 jets mass spectrum for anti-kt 4 and anti-kt 6 jets (array index is jet type from enum above)
    Histo1DGroupPtr  h_trijet_Mass[2];

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1326641);
}
