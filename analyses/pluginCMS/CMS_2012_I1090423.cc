// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2012_I1090423 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I1090423);

    void init() {

      FinalState fs;
      FastJets antikt(fs, JetAlg::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");

      book(_h_chi_dijet, {400., 600., 800., 1000., 1200., 1500., 1900., 2400., 3000., 7000.});
      for (auto& b : _h_chi_dijet->bins()) {
        book(b, _h_chi_dijet->numBins() - b.index() + 1, 1, 1);
      }
    }


    void analyze(const Event& event) {
      const Jets& jets = apply<JetFinder>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const double y0 = jets[0].rapidity();
      const double y1 = jets[1].rapidity();
      if (fabs(y0+y1)/2 > 1.11) vetoEvent;

      const double chi = exp(fabs(y0-y1));
      if (chi > 16) vetoEvent;

      const FourMomentum jj = jets[0].momentum() + jets[1].momentum();
       _h_chi_dijet->fill(jj.mass(), chi);
    }


    void finalize() {
      normalize(_h_chi_dijet);
    }


  private:

    Histo1DGroupPtr _h_chi_dijet;

  };



  RIVET_DECLARE_PLUGIN(CMS_2012_I1090423);

}
