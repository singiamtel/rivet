// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief EW W+W- production at 13 TeV
  class CMS_2022_I2080534 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2022_I2080534);

    /// @name Analysis methods
    /// @{
    void init() {

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9 && Cuts::pT > 100 * MeV);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, JetAlg::ANTIKT, 0.4);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 10 * GeV;
      LeptonFinder dressed_leps(bare_leps, photons, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      book(_h_xsec_exp, 1, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve dressed leptons, sorted by pT
      const Particles& leptons = apply<LeptonFinder>(event, "leptons").particlesByPt();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30 * GeV && Cuts::abseta < 4.7);

      // Retrieve Etmiss variable
      const double Etmiss = apply<MissingMomentum>(event, "MET").missingEt();

      // Remove jets with dressed leptons within dR < 0.4
      idiscardIfAnyDeltaRLess(jets, leptons, 0.4);
      if (jets.size() < 2)  vetoEvent;

      // Veto event if there is at least one b-jet
      for (const Jet& j : jets) {
        if (j.bTagged())  vetoEvent;
      }

      // Apply lepton cuts
      if (leptons.size() < 2)  vetoEvent;
      if (leptons[0].pT() < 25 * GeV || leptons[1].pT() < 13 * GeV)  vetoEvent;
      if (leptons.size() > 2 && leptons[2].pT() > 10 * GeV)          vetoEvent;
      if (leptons[0].pid() * leptons[1].pid() != -11 * 13 &&
          leptons[0].pid() * leptons[1].pid() != -11 * 11 &&
          leptons[0].pid() * leptons[1].pid() != -13 * 13)  vetoEvent;

      FourMomentum l1 = leptons[0].momentum();
      FourMomentum l2 = leptons[1].momentum();
      FourMomentum ll = (l1 + l2);
      const double mll = ll.mass();

      if (ll.pT() < 30*GeV)  vetoEvent;
      if (mll < 50*GeV)      vetoEvent;

      // Apply a missing-momentum cut
      if (Etmiss < 20*GeV)   vetoEvent;

      FourMomentum j1 = jets[0].momentum();
      FourMomentum j2 = jets[1].momentum();
      const double mjj = (j1 + j2).mass();
      const double detajj = abs(j1.eta() - j2.eta());

      // Apply VBS cuts
      if (mjj < 300 * GeV || detajj < 2.5)  vetoEvent;

      _h_xsec_exp->fill("Exclusive"s);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double norm = crossSection()  / femtobarn / sumOfWeights();

      scale(_h_xsec_exp, norm);
    }

    //@}

    /// @name Histograms
    //@{
    BinnedHistoPtr<string> _h_xsec_exp;
    //@}
  };

  RIVET_DECLARE_PLUGIN(CMS_2022_I2080534);

}  // namespace Rivet
