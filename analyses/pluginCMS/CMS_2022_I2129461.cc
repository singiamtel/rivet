// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Normalised differential cross sections for 13 TeV tW-channel single top-quark production
  class CMS_2022_I2129461 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2022_I2129461);


    /// Book histograms and initialise projections before the run
    void init() {

      // final state of all stable particles
      Cut particle_cut = (Cuts::abseta < 5.0) and (Cuts::pT > 0.*MeV);
      FinalState fs(particle_cut);

      // select charged leptons
      ChargedLeptons charged_leptons(fs);

      // select final state photons for dressed lepton clustering
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // select final state prompt charged leptons
      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);

      // select final state prompt photons
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);

      // Dressed leptons from selected prompt charged leptons and photons
      Cut lepton_cut = ((Cuts::abseta < 2.4) and 
                        (Cuts::pT > 20.*GeV) and 
                        (((Cuts::abspid == PID::ELECTRON) and ((Cuts::abseta < 1.4442) or (Cuts::abseta > 1.566))) or (Cuts::abspid == PID::MUON)));

      LeptonFinder dressed_leptons(
        prompt_photons, prompt_leptons, 0.1, lepton_cut);
      declare(dressed_leptons, "LeptonFinder");

      // Jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      declare(
        // excludes all neutrinos by default
        FastJets(fsForJets, JetAlg::ANTIKT, 0.4),
        "Jets");

      // pTmiss
      declare(MissingMomentum(fs), "MET");

      // book histograms
      book(_hist_norm_lep1_pt,            "d06-x01-y08");
      book(_hist_norm_lep1lep2jet1_pz,    "d07-x01-y08");
      book(_hist_norm_jet1_pt,            "d08-x01-y08");
      book(_hist_norm_lep1lep2jet1_m,     "d09-x01-y08");
      book(_hist_norm_lep1lep2_dphi,      "d10-x01-y08");
      book(_hist_norm_lep1lep2jet1met_mt, "d11-x01-y08");
    }


    /// @brief Perform the per-event analysis
    void analyze(const Event& event) {
      DressedLeptons dressedLeptons = apply<LeptonFinder>(
        event,
        "LeptonFinder"
      ).dressedLeptons();

      // Require at least two dressed leptons
      if (dressedLeptons.size() < 2) vetoEvent;

      // Require that the leading ones are an electron and a muon of opposite charge
      if (abs(dressedLeptons.at(0).pid() + dressedLeptons.at(1).pid()) != 2) vetoEvent;

      // Require that the leading lepton has at least 25 GeV of pT
      if (dressedLeptons.at(0).pt() <= 25.*GeV) vetoEvent;

      // Require that all identified lepton pairs have at least 20 GeV of invariant mass
      for (int iL = 0; iL < int(dressedLeptons.size()); iL++) {
        for (int jL = iL + 1; jL < int(dressedLeptons.size()); jL++) {
          if ((dressedLeptons.at(iL).momentum() + dressedLeptons.at(jL).momentum()).mass() <= 20.*GeV) vetoEvent;
        }
      }

      // Jet object ID
      Cut jet_cut((Cuts::abseta < 2.4) and (Cuts::pT > 20.*GeV));
      vector<Jet> jets = apply<FastJets>(
        event,
        "Jets"
      ).jets(jet_cut);

      // ignore jets that overlap with dressed leptons within dR < 0.4 and collect them in loose jets and jets categories
      Jets cleanedLooseJets;
      Jets cleanedJets;
      idiscardIfAnyDeltaRLess(jets, dressedLeptons, 0.4);
      for (const Jet& jet: jets) {
        if (jet.pt() > 30.*GeV) cleanedJets.push_back(jet);
        else                    cleanedLooseJets.push_back(jet);
      }

      if (cleanedJets.size() != 1)         vetoEvent; // select events with exactly one jet...
      if (not cleanedJets.at(0).bTagged()) vetoEvent; // ...that must be b-tagged...
      if (cleanedLooseJets.size() > 0)     vetoEvent; // ...and no loose jets

      // fill the histograms
      _hist_norm_lep1_pt->fill(min(max(dressedLeptons.at(0).pt(), 26.*GeV), 149.*GeV) / GeV);
      FourMomentum llj = dressedLeptons.at(0).momentum() + dressedLeptons.at(1).momentum() + cleanedJets.at(0).momentum();
      _hist_norm_lep1lep2jet1_pz->fill(min(max(abs((llj).pz()), 1.*GeV), 449.*GeV) / GeV);
      _hist_norm_jet1_pt->fill(min(max(cleanedJets.at(0).pt(), 31.*GeV), 149.*GeV) / GeV);
      _hist_norm_lep1lep2jet1_m->fill(min(max((llj).mass(), 51.*GeV), 399.*GeV) / GeV);
      _hist_norm_lep1lep2_dphi->fill( abs(deltaPhi(dressedLeptons.at(0).phi(), dressedLeptons.at(1).phi())) / PI );

      Vector3 tmpPtmiss3  = apply<MissingMomentum>(event, "MET").vectorPt();
      tmpPtmiss3.setZ(0.);
      FourMomentum ptmiss = FourMomentum(tmpPtmiss3.mod(), -tmpPtmiss3.x(), -tmpPtmiss3.y(), 0.);
      FourMomentum tmpM   = llj + ptmiss;
      _hist_norm_lep1lep2jet1met_mt->fill(min(max(sqrt(tmpM.E2() - tmpM.pz2()), 101.*GeV), 499.*GeV) / GeV);
    }

    /// @brief Normalise histograms after the run
    void finalize() {
      normalize(_hist_norm_lep1_pt);
      normalize(_hist_norm_lep1lep2jet1_pz);
      normalize(_hist_norm_jet1_pt);
      normalize(_hist_norm_lep1lep2jet1_m);
      normalize(_hist_norm_lep1lep2_dphi);
      normalize(_hist_norm_lep1lep2jet1met_mt);
    }

  private:
    // Declaration of histograms
    Histo1DPtr _hist_norm_lep1_pt;
    Histo1DPtr _hist_norm_lep1lep2jet1_pz;
    Histo1DPtr _hist_norm_jet1_pt;
    Histo1DPtr _hist_norm_lep1lep2jet1_m;
    Histo1DPtr _hist_norm_lep1lep2_dphi;
    Histo1DPtr _hist_norm_lep1lep2jet1met_mt;
  };

  RIVET_DECLARE_PLUGIN(CMS_2022_I2129461);
}
