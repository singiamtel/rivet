// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief VBS Zy production at 13 TeV
  class ATLAS_2023_I2663725 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2663725);

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      // Prompt photons
      const PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::pT > 25 * GeV && Cuts::abseta < 2.37);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState bareelectron_fs = Cuts::abspid == PID::ELECTRON;
      const PromptFinalState baremuon_fs = Cuts::abspid == PID::MUON;

      // Dressed leptons
      const FinalState allphoton_fs(Cuts::abspid == PID::PHOTON);
      const Cut leptoncut = Cuts::pT > 20 * GeV && Cuts::abseta < 2.47;
      // use *all* photons for dressing
      const LeptonFinder dressedelectron_fs(bareelectron_fs, allphoton_fs, 0.1, leptoncut);
      const LeptonFinder dressedmuon_fs(baremuon_fs, allphoton_fs, 0.1, leptoncut);

      declare(dressedelectron_fs, "Electrons");
      declare(dressedmuon_fs, "Muons");

      // FS excluding the leading photon
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photon_fs);
      vfs.addVetoOnThisFinalState(dressedmuon_fs);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "isolatedFS");

      VetoedFinalState hadrons(FinalState(Cuts::abseta < 4.4));
      hadrons.addVetoOnThisFinalState(dressedelectron_fs);
      hadrons.addVetoOnThisFinalState(dressedmuon_fs);

      // Jets
      FastJets jets(hadrons, JetAlg::ANTIKT, 0.4, JetMuons::ALL, JetInvisibles::DECAY);
      declare(jets, "jets");

      // Histograms
      size_t idx = 0;
      const vector<string> dists{"pTlep1", "pTy", "pTjet1", "pTll", "pTlly", "mjj", "dRap", "dPhi", "cent"};
      for (const string& reg : vector<string>{"", "_ext"}) {
        for (const string& obs : dists) {
          if (reg == "" && (obs == "pTll" || obs == "cent"))  continue;
          book(_h[obs+reg], ++idx, 1, 1);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event &event) {

      // Sort the dressed leptons by pt
      Particles electrons = apply<LeptonFinder>(event, "Electrons").particlesByPt();
      Particles muons = apply<LeptonFinder>(event, "Muons").particlesByPt();

      if (electrons.size() < 2 && muons.size() < 2) vetoEvent;

      Particles lep;
      if (electrons.size() >= 2) {
        lep.push_back(electrons[0]);
        lep.push_back(electrons[1]);
      }
      else if (muons.size() >= 2) {
        lep.push_back(muons[0]);
        lep.push_back(muons[1]);
      }
      if (lep.size() != 2)  vetoEvent;
      if (lep[0].pT() < 30*GeV) vetoEvent;

      const double mll = (lep[0].momentum() + lep[1].momentum()).mass();

      Particles photons = apply<PromptFinalState>(event, "Photons").particlesByPt();
      if (photons.empty())  vetoEvent;

      idiscardIfAnyDeltaRLess(photons, electrons, 0.4);
      idiscardIfAnyDeltaRLess(photons, muons, 0.4);

      Particles selectedPh;
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      for (const Particle& ph : photons) {
        // check photon isolation
        double coneEnergy(0.0);
        for (const Particle& p : fs) {
          if (deltaR(ph, p) < 0.2)  coneEnergy += p.Et(); // etcone20
        }
        if (coneEnergy / ph.pT() > 0.07)  continue;
        selectedPh += ph;
      }

      if (selectedPh.size() < 1)  vetoEvent;

      // Get jets
      const Cut jetscut = (Cuts::pT > 25 * GeV && Cuts::absrap < 4.4);
      Jets jetsMix = apply<FastJets>(event, "jets").jetsByPt(jetscut);

      idiscardIfAnyDeltaRLess(jetsMix, photons, 0.4);
      idiscardIfAnyDeltaRLess(jetsMix, electrons, 0.3);
      idiscardIfAnyDeltaRLess(jetsMix, muons, 0.3);

      // Jet multiplicity
      const size_t njets = jetsMix.size();
      if (njets < 2)  vetoEvent;

      if (jetsMix[1].pT() < 50*GeV) vetoEvent; // sorted by pT, must all > 50GeV

      if (mll < 40*GeV)   vetoEvent;

      const FourMomentum lly_system = lep[0].momentum() + lep[1].momentum() + selectedPh[0].momentum();
      const double mlly = lly_system.mass();
      if (mll + mlly <= 182 * GeV) vetoEvent;

      // other complex cuts
      const FourMomentum jj_system = jetsMix[0].momentum() + jetsMix[1].momentum();
      const double mjj = jj_system.mass();
      if (mjj < 150*GeV)  vetoEvent;

      const double Dyjj = fabs(deltaRap(jetsMix[0].momentum(), jetsMix[1].momentum()));
      if (Dyjj < 1)  vetoEvent;

      const double j1rap = jetsMix[0].rap();
      const double j2rap = jetsMix[1].rap();

      const double Centrality = fabs((lly_system.rap() - 0.5*(j1rap + j2rap)) / (j1rap - j2rap));
      if (Centrality > 5) vetoEvent;

      size_t njGap = 0;
      for (size_t i=2; i < njets; ++i) {
        const double jrap = jetsMix[i].rap();
        if ((jrap < j1rap && jrap > j2rap) || (jrap < j2rap && jrap > j1rap))  ++njGap;
      }
      if (njGap > 0)  vetoEvent;

      // now start to log hist
      // first we select the extended SR
      if (Centrality > 0.4) vetoEvent;

      _h["pTjet1_ext"]->fill(jetsMix[0].pT() / GeV);
      _h["pTy_ext"]->fill(selectedPh[0].pT() / GeV);
      _h["pTlep1_ext"]->fill(lep[0].pT() / GeV);
      _h["pTll_ext"]->fill((lep[0].momentum() + lep[1].momentum()).pT() / GeV);
      _h["pTlly_ext"]->fill(lly_system.pT() / GeV);

      _h["cent_ext"]->fill(Centrality);
      _h["dRap_ext"]->fill(Dyjj);
      _h["mjj_ext"]->fill(mjj);
      _h["dPhi_ext"]->fill(fabs(deltaPhi(lly_system, jj_system)));

      // deal with EWK SR500
      if (mjj < 450 * GeV)  vetoEvent; // special since mjj in SR includes the underflow

      _h["mjj"]->fill(mjj);

      if (mjj < 500*GeV)  vetoEvent;

      _h["pTjet1"]->fill(jetsMix[0].pT() / GeV);
      _h["pTy"]->fill(selectedPh[0].pT() / GeV);
      _h["pTlep1"]->fill(lep[0].pT() / GeV);
      _h["pTlly"]->fill(lly_system.pT() / GeV);
      _h["dRap"]->fill(Dyjj);
      _h["dPhi"]->fill(fabs(deltaPhi(lly_system, jj_system)));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection() / femtobarn / sumOfWeights());
    }

  private:

    /// Histograms
    map<string, Histo1DPtr> _h;

  };

  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2663725);

}
