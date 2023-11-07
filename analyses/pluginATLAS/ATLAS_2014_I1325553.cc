// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Jet mass as a function of ystar
  class ATLAS_2014_I1325553 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1325553);

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

      const vector<double> ybins{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      size_t ptDsOffset = 0;
      for (size_t alg = 0; alg < 2; ++alg) {
        book(_pt[alg], ybins);
        for (auto& b : _pt[alg]->bins()) {
          book(b, 1+ptDsOffset, 1, 1);
          ++ptDsOffset;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jetAr[2];
      jetAr[AKT4] = apply<FastJets>(event, "AntiKT04").jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 3.0);
      jetAr[AKT6] = apply<FastJets>(event, "AntiKT06").jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 3.0);

      // Loop over jet "radii" used in analysis
      for (size_t alg = 0; alg < 2; ++alg) {

        // fill the 1D pt histograms with all the jets passing the cuts
        for (const Jet& jet : jetAr[alg]) {
          const double absrap = jet.absrap();
          if (absrap < 3.0) {
            const double pt = jet.pT();
            if (pt/GeV > 100*GeV) {
              _pt[alg]->fill(absrap, pt/GeV);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Print summary info
      const double xs_norm_factor = 0.5*crossSection()/picobarn/sumOfWeights();
      scale(_pt, xs_norm_factor);
      divByGroupWidth(_pt);
    }

    /// @}


  private:

    // Data members like post-cuts event weight counters go here
    enum Alg { AKT4=0, AKT6=1 };

    /// The inclusive jet spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    Histo1DGroupPtr _pt[2];

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1325553);

}
