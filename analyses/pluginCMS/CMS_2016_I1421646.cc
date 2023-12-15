// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS azimuthal decorrelations at 8 TeV
  class CMS_2016_I1421646 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_I1421646);


    /// Book projections and histograms
    void init() {

      FastJets akt(FinalState(), JetAlg::ANTIKT, 0.7);
      declare(akt, "antikT");

      book(_h_deltaPhi, {200., 300., 400., 500., 700., 900., 1100., 4000.});
      for (auto& b : _h_deltaPhi->bins()) {
        book(b, b.index(), 1, 1);
      }
    }


    /// Per-event analysis
    void analyze(const Event & event) {

      const Jets& jets = apply<JetFinder>(event, "antikT").jetsByPt(Cuts::absrap < 5.0 && Cuts::pT > 100*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].pT() < 200*GeV) vetoEvent;
      if (jets[0].absrap() > 2.5 || jets[1].absrap() > 2.5) vetoEvent;

      const double dphi = deltaPhi(jets[0].phi(), jets[1].phi());
      _h_deltaPhi->fill(jets[0].pT(), dphi);
    }


    /// Scale histograms
    void finalize() {
      normalize(_h_deltaPhi);
    }


  private:

    Histo1DGroupPtr _h_deltaPhi;

  };


  RIVET_DECLARE_PLUGIN(CMS_2016_I1421646);

}
