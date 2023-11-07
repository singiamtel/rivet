// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2015_I1327224 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1327224);


    void init() {

      FinalState fs;
      FastJets antikt(fs, JetAlg::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");

      book(_h_chi_dijet, {1900., 2400., 3000., 3600., 4200., 8000.});
      for (auto& b : _h_chi_dijet->bins()) {
        book(b, _h_chi_dijet->numBins() - b.index() + 1, 1, 1);
      }
    }


    void analyze(const Event& event) {

      const Jets& jets = apply<JetFinder>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      const double y0 = j0.rapidity();
      const double y1 = j1.rapidity();
      if (fabs(y0 + y1) / 2. > 1.11) vetoEvent;

      const double mjj = FourMomentum(j0 + j1).mass();
      if (mjj/GeV <1900) vetoEvent;

      const double chi = exp(fabs(y0 - y1));
      if (chi >= 16.)  vetoEvent;

      // Fill the histogram
      _h_chi_dijet->fill(mjj/GeV, chi);
    }

    void finalize() {
      normalize(_h_chi_dijet);
    }

  private:

    Histo1DGroupPtr _h_chi_dijet;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2015_I1327224);
}
