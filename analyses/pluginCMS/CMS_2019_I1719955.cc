
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// CMS azimuthal decorrelations in back-to-back dijet events at 13 TeV
  class CMS_2019_I1719955 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2019_I1719955);


    /// Book projections and histograms
    void init() {
      const FinalState fs;
      declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "ANTIKT");

      book(_h_deltaPhi_2J, {200., 300., 400., 500., 600., 700., 800., 1000., 1200., 4000.});
      book(_h_deltaPhi_3J, {200., 300., 400., 500., 600., 700., 800., 1000., 1200., 4000.});
      for (size_t i=1; i<_h_deltaPhi_2J->numBins()+1; ++i) {
        book(_h_deltaPhi_2J->bin(i), i, 1, 1);
        book(_h_deltaPhi_3J->bin(i), i+9, 1, 1);
      }
    }

    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<JetFinder>(event, "ANTIKT").jetsByPt(Cuts::absrap < 5. && Cuts::pT > 100*GeV);
      const Jets& lowjets = apply<JetFinder>(event, "ANTIKT").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 30*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].absrap() > 2.5 || jets[1].absrap() > 2.5) vetoEvent;

      const double dphi = 180./M_PI*deltaPhi(jets[0].phi(), jets[1].phi());
      _h_deltaPhi_2J->fill(jets[0].pT(), dphi);
      if (lowjets.size() > 2) _h_deltaPhi_3J->fill(jets[0].pT(), dphi);
    }

    /// Scale histograms
    void finalize() {
      int region_ptmax_2J = 0;
      double norm_finalize[9];
      for (auto& histo_2J : _h_deltaPhi_2J->bins()) {
        norm_finalize[region_ptmax_2J] = histo_2J->integral();
        if (norm_finalize[region_ptmax_2J] != 0) scale(histo_2J, 1.0/norm_finalize[region_ptmax_2J]);
        region_ptmax_2J++;
      }
      int region_ptmax_3J = 0;
      for (auto& histo_3J : _h_deltaPhi_3J->bins()) {
        if (norm_finalize[region_ptmax_3J] != 0) scale(histo_3J, 1.0/norm_finalize[region_ptmax_3J]);
        region_ptmax_3J++;
      }
    }


  private:

    Histo1DGroupPtr _h_deltaPhi_2J, _h_deltaPhi_3J;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2019_I1719955);

}
