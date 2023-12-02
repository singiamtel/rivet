// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


  /// @brief Colour reconnection in ttbar dilepton events
  class ATLAS_2022_I2152933 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2152933);


    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

    	// Photons
    	PromptFinalState photons(Cuts::abspid == PID::PHOTON);

    	// Muons
    	Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
    	PromptFinalState bare_mu(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
    	LeptonFinder all_dressed_mu(bare_mu, photons, 0.1, lepton_cuts);
    	declare(all_dressed_mu, "muons");

    	// Electrons
    	PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
    	LeptonFinder all_dressed_el(bare_el, photons, 0.1, lepton_cuts);
    	declare(all_dressed_el, "electrons");

    	//Jet forming
    	const InvisibleFinalState neutrinos(OnlyPrompt::YES, TauDecaysAs::PROMPT);

    	VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0));
    	vfs.addVetoOnThisFinalState(all_dressed_el);
    	vfs.addVetoOnThisFinalState(all_dressed_mu);
    	vfs.addVetoOnThisFinalState(neutrinos);
    	FastJets jetfs(vfs, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::ALL);
    	declare(jetfs, "jets");

      // FinalState charged particles
      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");

      declare(MissingMomentum(), "ETmiss");

      // Book histograms
      book(_h["nch"], 7,1,1);
      book(_h["sumPt"], 8,1,1);

      book(_h_multi, { 0, 19.5, 39.5, 59.5, 79.5, 101 });
      for (auto& b : _h_multi->bins()) {
        book(b, 8+b.index(), 1, 1);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve dressed leptons, sorted by pT
      DressedLeptons electrons = apply<LeptonFinder>(event, "electrons").dressedLeptons();
      DressedLeptons muons = apply<LeptonFinder>(event, "muons").dressedLeptons();
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      // Retrieve charge final state particles, sorted by pT
      const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();

      // OVERLAP REMOVAL
      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
      Jets bjets = select(jets, [](const Jet& jet) {
        return  jet.bTagged(Cuts::pT > 5*GeV);
      });

      // Veto event if there are not exactly one electron and one muon
      if ( electrons.size() != 1 || muons.size() != 1)  vetoEvent;

      // Veto event if the selected electron and muon are not OS
      if (electrons[0].charge() == muons[0].charge())  vetoEvent;

      // Veto event if there are no 2 jets or 3 jets
      if (jets.size() < 2 || jets.size() > 3)  vetoEvent;

      // Veto event if the number of b-jets not exactly 2
      if (bjets.size() != 2)  vetoEvent;

      // Veto event with dilepton mass < 15 GeV
      FourMomentum ll = electrons[0].momentum() + muons[0].momentum();
      if (ll.mass() <= 15*GeV ) vetoEvent;

      // Remove tracks inside jets and leptons: within 0.4 delatR from a jet and 0.01 from a lepton

      Particles myTracks;
      for (const Particle& p : tracks) {
        bool isJetAssoc = false;
        for (unsigned int i = 0; i < jets.size(); ++i) {
          if(deltaR(jets[i], p, PSEUDORAPIDITY) < 0.4) {
          	isJetAssoc= true;
          	break;
          }
        }
        if (isJetAssoc) continue;
        if(deltaR(muons[0], p, PSEUDORAPIDITY) < 0.01) continue;
        if(deltaR(electrons[0], p, PSEUDORAPIDITY) < 0.01) continue;
        if(p.pid() == 0 || p.isStable() != 1 || p.charge() == 0) continue;
        myTracks.push_back(p);
      }

      // Veto event if the number of selected charged particles higher than 100
      if (myTracks.size() > 100)  vetoEvent;

      // Fill nch, sumpt, and sumpt vs. nch
      const double nch = min(myTracks.size(), (size_t)99); // put overflow in the last bin
      _h["nch"]->fill(nch);

      const double sumPt = min(sum(myTracks, Kin::pT, 0.0), 119*GeV);
      _h["sumPt"]->fill(sumPt/GeV);

      const double sumPt2d = min(sumPt, 99*GeV);
      _h_multi->fill(nch, sumPt2d/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize to unity
      normalize(_h);
      normalizeGroup(_h_multi, 1.0, false);
      divByGroupWidth(_h_multi);
    }


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;
    Histo1DGroupPtr _h_multi;

    /// @}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2152933);

}
