// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// W + jets jet multiplicities and pT
  class ATLAS_2010_S8919674 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2010_S8919674);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare("MET", MissingMomentum());

      // Find electrons and muons
      const Cut cut_e = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47)) && Cuts::pT > 20*GeV;
      LeptonFinder ef(0.1, cut_e && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");
      const Cut cut_m = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      LeptonFinder mf(0.1, cut_m && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      // Input for the jets: no neutrinos, no muons, and no electron which passed the electron cuts
      VetoedFinalState veto;
      veto.addVetoOnThisFinalState(ef);
      veto.addVetoOnThisFinalState(mf);
      FastJets jets(veto, JetAlg::ANTIKT, 0.4, JetMuons::NONE);
      declare(jets, "Jets");

      /// Book histograms
      book(_h_el_njet_inclusive,1,1,1);
      book(_h_mu_njet_inclusive,2,1,1);
      book(_h_el_pT_jet1,5,1,1);
      book(_h_mu_pT_jet1,6,1,1);
      book(_h_el_pT_jet2,7,1,1);
      book(_h_mu_pT_jet2,8,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _h_mu_njet_inclusive->xEdges();

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV);

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.Et() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 1000*GeV);

      // Require two valid W candidates
      if (iefound < 0 && imfound < 0) vetoEvent; //< no W's
      if (iefound >= 0 && imfound >= 0) vetoEvent; //< multi-W

      // Histogramming
      if (iefound >= 0) {
        const Particle& e = es[iefound].constituents()[0];
        if (mT(pmiss, e) > 40*GeV) {
          const Jets js = select(jets, [&](const Jet& j) { return j.abseta() < 2.8 && deltaR(e, j) > 0.5; });
          _h_el_njet_inclusive->fill(_edges[0]);
          if (js.size() >= 1) {
            _h_el_njet_inclusive->fill(_edges[1]);
            _h_el_pT_jet1->fill(js[0].pT());
          }
          if (js.size() >= 2) {
            _h_el_njet_inclusive->fill(_edges[2]);
            _h_el_pT_jet2->fill(js[1].pT());
          }
          if (js.size() >= 3) {
            _h_el_njet_inclusive->fill(_edges[3]);
          }
        }
      }

      if (imfound >= 0) {
        const Particle& mu = mus[imfound];
        if (mT(pmiss, mu) > 40*GeV) {
          const Jets js = select(jets, [&](const Jet& j) { return j.abseta() < 2.8 && deltaR(mu, j) > 0.5; });
          _h_mu_njet_inclusive->fill(_edges[0]);
          if (js.size() >= 1) {
            _h_mu_njet_inclusive->fill(_edges[1]);
            _h_mu_pT_jet1->fill(js[0].pT());
          }
          if (js.size() >= 2) {
            _h_mu_njet_inclusive->fill(_edges[2]);
            _h_mu_pT_jet2->fill(js[1].pT());
          }
          if (js.size() >= 3) {
            _h_mu_njet_inclusive->fill(_edges[3]);
          }
          if (js.size() >= 4) {
            _h_mu_njet_inclusive->fill(_edges[4]);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double normfac = crossSection()/sumOfWeights();
      scale(_h_el_njet_inclusive, normfac/nanobarn);
      scale(_h_mu_njet_inclusive, normfac/nanobarn);
      scale(_h_el_pT_jet1, normfac);
      scale(_h_mu_pT_jet1, normfac);
      scale(_h_el_pT_jet2, normfac);
      scale(_h_mu_pT_jet2, normfac);
    }

    ///@}


  private:

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_el_njet_inclusive;
    BinnedHistoPtr<string> _h_mu_njet_inclusive;
    Histo1DPtr _h_el_pT_jet1;
    Histo1DPtr _h_mu_pT_jet1;
    Histo1DPtr _h_el_pT_jet2;
    Histo1DPtr _h_mu_pT_jet2;
    vector<string> _edges;
    /// @}

  };


  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2010_S8919674, ATLAS_2010_I882534);

}
