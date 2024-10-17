// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief ttbar + gamma at 13 TeV
  class ATLAS_2024_I2768921 : public Analysis {
  public:


    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2024_I2768921);


    // Book histograms and initialise projections before the run
    void init() {

      // Set running mode (default is both channels)
      _mode = 3;
      if ( getOption("LMODE") == "SINGLE" )   _mode = 1;
      if ( getOption("LMODE") == "DILEPTON" ) _mode = 2;
      if ( getOption("LMODE") == "ALL" )      _mode = 3;

      // Cuts
      Cut eta_full = Cuts::abseta < 5.0;
      Cut lepCuts  = Cuts::abseta < 2.5 && Cuts::pT > 25.0*GeV;
      Cut lepCutsVeto = Cuts::abseta < 2.5 && Cuts::pT > 7.0*GeV;

      // Charged particles for signal photon isolation
      declare(ChargedFinalState(), "CFS");

      // Signal photons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON && Cuts::pT > 15*GeV && Cuts::abseta < 2.5);
      declare(photons, "Photons");

      // Invisibles
      InvisibleFinalState prompt_invis(OnlyPrompt::YES, TauDecaysAs::PROMPT);

      // Leptons
      PromptFinalState leptons(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);

      // Dress the leptons
      FinalState dressPhotons(Cuts::abspid == PID::PHOTON);
      LeptonFinder dressedLeptons(leptons, dressPhotons, 0.1, lepCuts);
      declare(dressedLeptons, "DressedLeptons");
      // Full dressed leptons
      LeptonFinder fulldressedleptons(leptons, dressPhotons, 0.1, eta_full);

      // Add new lepton collection for veto of 2nd lepton pt > 7 GeV
      LeptonFinder vetodressedLeptons(leptons, dressPhotons, 0.1, lepCutsVeto);
      declare(vetodressedLeptons, "VetoDressedLeptons");

      // Jet alg input
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(fulldressedleptons);
      vfs.addVetoOnThisFinalState(prompt_invis);

      // Jet clustering
      FastJets jets(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
      declare(jets, "Jets");

      dualbook("sldl_ph_pt",     18);
      dualbook("sldl_ph_abseta", 20);

      // Book histograms
      if (_mode == 1 || _mode == 3) { // Single lepton channel
        dualbook("sl_ph_bjet_dR",  6);
        dualbook("sl_ph_l_dR",     8);
        dualbook("sl_l_jet_dR",   10);
        dualbook("sl_ph_pt",      44);
        dualbook("sl_ph_abseta",  46);
        dualbook("sl_jet1_pt",    48);
      }

      if (_mode == 2 || _mode == 3) { // Dilepton channel
        dualbook("dl_ph_bjet_dR", 12);
        dualbook("dl_ph_l_dR",    14);
        dualbook("dl_l_jet_dR",   16);
        dualbook("dl_ph_pt",      50);
        dualbook("dl_ph_abseta",  52);
        dualbook("dl_jet1_pt",    54);
      }
    }


    void analyze(const Event& event) {

      FourMomentum lepton(0.,0.,0.,0.);

      // Fetch objects
      Particles photons          = apply<PromptFinalState>(event, "Photons").particlesByPt();
      DressedLeptons leptons     = apply<LeptonFinder>(event, "DressedLeptons").dressedLeptons();
      DressedLeptons vetoleptons = apply<LeptonFinder>(event, "VetoDressedLeptons").dressedLeptons();
      Jets jets                  = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);


      // Immediate veto on events without photon
      if (photons.empty())  vetoEvent;

      // Overlap removal of leptons near jets
      idiscardIfAnyDeltaRLess(leptons, jets, 0.4);
      idiscardIfAnyDeltaRLess(vetoleptons, jets, 0.4);

      // Photon isolation with overlap removal of jets near isolated photon
      Particles iso_photons;
      const ChargedFinalState& chargedFS = apply<ChargedFinalState>(event, "CFS");
      for (const Particle& ph : photons) {
        const double conePt = sum(chargedFS.particles(deltaRLess(ph, 0.2)), Kin::pT, 0.0);
        if (conePt/ph.pT() < 0.05) {
          iso_photons += ph;
          idiscard(jets, deltaRLess(ph, 0.4));
        }
      }

      // Veto on events without one good photon
      if (iso_photons.size() != 1) vetoEvent;

      // Require at least one bjet
      Jets bjets;
      for (const Jet& jet : jets) {
       if (jet.bTagged(Cuts::pT > 5*GeV)) bjets += jet;
      }
      if (bjets.empty()) vetoEvent;


      const Particle& photon = iso_photons[0];
      // Apply cuts to isolated photon, after requiring exactly 1 photon
      if (photon.pT() <= 20*GeV || photon.abseta() >= 2.37) vetoEvent;

      // Veto event if photon too close to a lepton
      if ( any(leptons, deltaRLess(photon, 0.4)) )  vetoEvent;


      // Fill histos in single lepton channel, last bin contains overflow events
      if (_mode == 1 || _mode == 3) {
        if (jets.size() >= 4 && leptons.size() == 1 && vetoleptons.size() < 2) {

          double mindR_ph_l = 1e3;
          for (const DressedLepton& lep : leptons) {
            double dr = deltaR(photon, lep);
            if (dr <  mindR_ph_l)  mindR_ph_l = dr;
          }
          if (mindR_ph_l >= 3.4)  mindR_ph_l = 3.;

          // dRmin(photon,bjet)
          double mindR_ph_bjet = 1e3;
          for (const Jet& bjet : bjets) {
            double dr = deltaR(bjet, photon);
            if (dr < mindR_ph_bjet) mindR_ph_bjet = dr;
          }
          if (mindR_ph_bjet >= 3.4) mindR_ph_bjet = 3.;

          //dRmin(lepton,jet)
          double mindR_l_jet = 1e3;
          for (const DressedLepton& lep : leptons) {
            for (const Jet& jet : jets) {
              double dr = deltaR(lep, jet);
              if (dr < mindR_l_jet) mindR_l_jet = dr;
            }
          }
          if (mindR_l_jet >= 3.4) mindR_l_jet = 3.;

          double ph_pt = photon.pT();
          if (ph_pt >= 500) ph_pt = 450.;

          double jet1_pt = jets[0].pT();
          if (jet1_pt >= 450.) jet1_pt = 405;

          dualfill("sl_ph_pt",       ph_pt/GeV);
          dualfill("sl_ph_abseta",   photon.abseta());
          dualfill("sl_ph_l_dR",     mindR_ph_l);
          dualfill("sl_ph_bjet_dR",  mindR_ph_bjet);
          dualfill("sl_l_jet_dR",    mindR_l_jet);
          dualfill("sl_jet1_pt",     jet1_pt/GeV);
          dualfill("sldl_ph_pt",     ph_pt/GeV);
          dualfill("sldl_ph_abseta", photon.abseta());

        }
      }


      // Fill histos in dilepton channel, last bin contains overflow events
      if ( _mode == 2 || _mode == 3 ) {
        if ( jets.size() >= 2 && leptons.size() == 2 && vetoleptons.size() < 3) {

          // dRmin(photon, lepton)
          double mindR_ph_l = 1e3;
          for (const DressedLepton& lep : leptons) {
            double dr = deltaR(lep, photon);
            if (dr < mindR_ph_l)  mindR_ph_l = dr;
          }
          if (mindR_ph_l >= 3.4)  mindR_ph_l = 3.;

          // dRmin(photon,bjet)
          double mindR_ph_bjet = 1e3;
          for (const Jet& bjet : bjets) {
            double dr = deltaR(bjet, photon);
            if (dr < mindR_ph_bjet ) mindR_ph_bjet = dr;
          }
          if (mindR_ph_bjet >= 3.4) mindR_ph_bjet = 3.;

          //dRmin(lepton,jet)
          double mindR_l_jet = 999.0;
          for (const DressedLepton& lep :leptons) {
            for (const Jet& jet : jets) {
             double dr = deltaR(lep, jet);
             if (dr < mindR_l_jet)  mindR_l_jet = dr;
            }
          }
          if (mindR_l_jet >= 3.4) mindR_l_jet = 3.;

          double ph_pt = photon.pT();
          if (ph_pt >= 500.) ph_pt = 450.;

          double jet1_pt = jets[0].pT();
          if (jet1_pt >= 450.) jet1_pt = 405.;

          dualfill("dl_ph_pt",       ph_pt/GeV);
          dualfill("dl_ph_abseta",   photon.abseta());
          dualfill("dl_ph_l_dR",     mindR_ph_l);
          dualfill("dl_ph_bjet_dR",  mindR_ph_bjet);
          dualfill("dl_l_jet_dR",    mindR_l_jet);
          dualfill("dl_jet1_pt",     jet1_pt/GeV);
          dualfill("sldl_ph_pt",     ph_pt/GeV);
          dualfill("sldl_ph_abseta", photon.abseta());

        }
      }
    }



    void finalize() {

      // Normalise histogram to absolute fiducial cross section or to unity
      scale(_h, crossSection() / femtobarn / sumOfWeights());
      normalize(_n, 1.0, true);

    }


    void dualbook(const string& name, unsigned int id) {
      book(_h["abs_"+name],    id, 1, 1);
      book(_n["norm_"+name], 1+id, 1, 1);
    }

    void dualfill(const string& name, const double value) {
      _h["abs_"+name]->fill(value);
      _n["norm_"+name]->fill(value);
    }

  private:

    /// @name Histograms
    map<string, Histo1DPtr> _h, _n;

    size_t _mode;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2024_I2768921);
}
