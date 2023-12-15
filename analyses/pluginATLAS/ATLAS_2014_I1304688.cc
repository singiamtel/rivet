// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ATLAS 7 TeV jets in ttbar events analysis
  ///
  /// @author W. H. Bell <W.Bell@cern.ch>
  /// @author A. Grohsjean <alexander.grohsjean@desy.de>
  class ATLAS_2014_I1304688 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1304688);

    void init() {
      // Eta ranges
      /// @todo 1 MeV? Really?
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;
      Cut eta_lep = Cuts::abseta < 2.5;

      // Get photons to dress leptons
      FinalState photons(eta_full && Cuts::abspid == PID::PHOTON);

      // Projection to find the electrons
      PromptFinalState electrons(eta_full && Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      LeptonFinder dressedelectrons(electrons, photons, 0.1, eta_lep && Cuts::pT > 25*GeV);
      declare(dressedelectrons, "dressedelectrons");
      LeptonFinder vetodressedelectrons(electrons, photons, 0.1, eta_lep && Cuts::pT >= 15*GeV);
      declare(vetodressedelectrons, "vetodressedelectrons");
      LeptonFinder ewdressedelectrons(electrons, photons, 0.1, eta_full);
      declare(ewdressedelectrons, "ewdressedelectrons");

      // Projection to find the muons
      PromptFinalState muons(eta_full && Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      LeptonFinder dressedmuons(muons, photons, 0.1, eta_lep && Cuts::pT >= 25*GeV);
      declare(dressedmuons, "dressedmuons");
      LeptonFinder vetodressedmuons(muons, photons, 0.1, eta_lep && Cuts::pT >= 15*GeV);
      declare(vetodressedmuons, "vetodressedmuons");
      LeptonFinder ewdressedmuons(muons, photons, 0.1, eta_full);
      declare(ewdressedmuons, "ewdressedmuons");

      // Projection to find neutrinos and produce MET
      InvisibleFinalState neutrinos(OnlyPrompt::YES, TauDecaysAs::PROMPT);
      declare(neutrinos, "neutrinos");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, JetAlg::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      // Book histograms
      for (size_t i = 0; i < pTcuts.size(); ++i) {
        const string name = "mult_"+std::to_string(i);
        book(_s[name], i+1, 1, 1);
      }
      for (size_t i = 0; i < njets; ++i) {
        const string name = "jet_"+std::to_string(i);
        book(_h[name], i+5, 1, 1);
      }
    }


    void analyze(const Event& event) {

      if (_edges.empty()) {
        _edges.resize(_s.size());
        for (size_t i = 0; i < _edges.size(); ++i) {
          const string name = "mult_"+std::to_string(i);
          _edges[i] = _s[name]->xEdges();
        }
      }

      // Get the selected objects, using the projections.
      const DressedLeptons dressedelectrons = sortByPt(apply<LeptonFinder>(event, "dressedelectrons").dressedLeptons());
      const DressedLeptons vetodressedelectrons = apply<LeptonFinder>(event, "vetodressedelectrons").dressedLeptons();

      const DressedLeptons dressedmuons = sortByPt(apply<LeptonFinder>(event, "dressedmuons").dressedLeptons());
      const DressedLeptons vetodressedmuons = apply<LeptonFinder>(event, "vetodressedmuons").dressedLeptons();

      if (dressedelectrons.empty() && dressedmuons.empty())  vetoEvent;
      if (dressedelectrons.size()) {
        if (vetodressedelectrons.size() > 1 || vetodressedmuons.size())  vetoEvent;
      }
      if (dressedmuons.size()) {
        if (vetodressedmuons.size() > 1 || vetodressedelectrons.size())  vetoEvent;
      }

      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      if (jets.size() < 3)  vetoEvent;


      // Calculate the missing ET, using the prompt neutrinos only (really?)
      /// @todo Why not use MissingMomentum?
      FourMomentum pmet;
      for (const Particle& p : apply<InvisibleFinalState>(event, "neutrinos").particlesByPt())  pmet += p.momentum();
      const double met_et = pmet.pT();
      const double met_phi = pmet.phi();
      if (met_et <= 30*GeV)  vetoEvent;
      if (dressedelectrons.size()) {
        if (transMass(dressedelectrons[0].pT(), dressedelectrons[0].phi(), met_et, met_phi) <= 35*GeV)  vetoEvent;
      }
      if (dressedmuons.size()) {
        if (transMass(dressedmuons[0].pT(), dressedmuons[0].phi(), met_et, met_phi) <= 35*GeV)  vetoEvent;
      }

      // Check overlap of jets/leptons.
      size_t jet_ntag = 0;
      bool overlap = false;
      for (size_t j1 = 0; j1 < jets.size(); ++j1) {
        const Jet& jet = jets[j1];
        // If dR(el,jet) < 0.4 skip the event
        for (const DressedLepton& el : dressedelectrons) {
          if (deltaR(jet, el) < 0.4)  overlap = true;
        }
        // If dR(mu,jet) < 0.4 skip the event
        for (const DressedLepton& mu : dressedmuons) {
          if (deltaR(jet, mu) < 0.4)  overlap = true;
        }
        // If dR(jet,jet) < 0.5 skip the event
        for (size_t j2 = j1+1; j2 < jets.size(); ++j2) {
          const Jet& jet2 = jets[j2];
          if (deltaR(jet, jet2) < 0.5)  overlap = true;
        }
        // Count the number of b-tags
        if (jet.bTags().size())  ++jet_ntag;
      }

      // Remove events with object overlap
      if (overlap || !jet_ntag)  vetoEvent;

      // Count the jet multiplicity for 25, 40, 60 and 80GeV
      vector<size_t> jet_n; jet_n.resize(pTcuts.size());
      for (size_t i = 0; i < pTcuts.size(); ++i) {
        jet_n[i] = countJets(jets, i);
        const string name = "mult_" + std::to_string(i);
        const string& edge = _edges[i][jet_n[i]-3];
        _s[name]->fill(edge);
      }

      // Check if the additional pT threshold requirements are passed
      const bool pass_jetPt = jets.size() > 1 && jets[0].pT() > 50*GeV && jets[1].pT() > 35*GeV;
      if (!pass_jetPt)  vetoEvent;


      // Fill histograms
      for (size_t i = 0; i < njets; ++i) {
        const string name = "jet_" + std::to_string(i);
        if (i > 1 && jet_n[0] <= i)  continue;
        _h[name]->fill(jets[i].pT()/GeV);
      }
    }


    void finalize() {
      // Normalize to cross-section x 0.5 to average lepton flavours
      const double norm = 0.5*crossSection()/sumOfWeights();
      scale(_h, norm);
      scale(_s, norm);
    }



  private:


    /// @name Physics object helper functions
    /// @{

    double transMass(double ptLep, double phiLep, double met, double phiMet) const {
      return sqrt(2.0*ptLep*met*(1 - cos(phiLep-phiMet)));
    }

    size_t countJets(const Jets& jets, size_t thresh) const {
      size_t jet_n = 0;
      for (const Jet& jet : jets) {
        if (jet.pT() > pTcuts[thresh])  ++jet_n;
      }
      const vector<size_t> ncutoff{8,7,6,5};
      return min(jet_n, ncutoff[thresh]);
    }

    /// @}


  private:

    /// @name Objects that are used by the event selection decisions
    /// @{

    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<string>> _s;
    vector<vector<string>> _edges;
    const size_t njets = 5;
    const vector<int> pTcuts{25,40,60,80};
    /// @}

  };



  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1304688);

}
