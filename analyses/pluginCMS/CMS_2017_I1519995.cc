// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Search for new physics with dijet angular distributions at 13 TeV
  class CMS_2017_I1519995 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1519995);


    /// Book projections and histograms
    void init() {

      FastJets antikt(FinalState(), FastJets::ANTIKT, 0.4);
      declare(antikt, "ANTIKT");

      book(_h_chi_dijet, {1900., 2400., 3600., 4200., 4800., 8000.});
      for (auto& b : _h_chi_dijet->bins()) {
        book(b, _h_chi_dijet->numBins() - b.index() + 1, 1, 1);
      }
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const FourMomentum j0(jets[0].mom()), j1(jets[1].mom());
      if (fabs(j0.rap()+j1.rap())/2 > 1.11) vetoEvent;

      const double mjj = (j0+j1).mass();
      const double chi = exp(fabs(j0.rap()-j1.rap()));
      if (chi < 16) _h_chi_dijet->fill(mjj/GeV, chi);
    }


    /// Normalize histograms
    void finalize() {
      normalize(_h_chi_dijet);
    }


   Histo1DGroupPtr _h_chi_dijet;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2017_I1519995);

}
