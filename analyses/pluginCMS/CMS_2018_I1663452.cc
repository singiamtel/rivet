// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2018_I1663452 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1663452);

    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.4);
      declare(antikt, "ANTIKT");
      book(_h_chi_dijet, {2400., 3000., 4200., 4800., 5400., 6000.});
      for (auto& b : _h_chi_dijet->bins()) {
        book(b, _h_chi_dijet->numBins() - b.index() + 1, 1, 1);
      }
    }

    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;
      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0+y1)/2. > 1.11) vetoEvent;
      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet->fill(mjj, chi);
    }


    void finalize() {
      normalize(_h_chi_dijet);
    }


  private:

    Histo1DGroupPtr _h_chi_dijet;

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2018_I1663452);

}
