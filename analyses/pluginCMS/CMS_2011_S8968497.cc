// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Dijet angular distributions and search for quark compositeness at 7 TeV
  class CMS_2011_S8968497 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_S8968497);


    void init() {

      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");

      book(_h_chi_dijet, {250., 350., 500., 650., 850., 1100., 1400., 1800., 2200., 7000.});
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



  RIVET_DECLARE_ALIASED_PLUGIN(CMS_2011_S8968497, CMS_2011_I889175);

}
