// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief WWgamma production in emu final state at 13 TeV
  class CMS_2023_I2709669 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2023_I2709669);

    /// @name Analysis methods
    /// @{
    void init() {
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance

      const FinalState fs(Cuts::abseta < 4.9);

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_electrons(Cuts::abspid == PID::ELECTRON, TauDecaysAs::PROMPT);
      PromptFinalState bare_muons(Cuts::abspid == PID::MUON, TauDecaysAs::PROMPT);
      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons

      Cut muon_cuts = Cuts::abseta < 2.4 && Cuts::pT > 20 * GeV;
      Cut electron_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25 * GeV;
      LeptonFinder dressed_electrons(bare_electrons, photons, 0.1, electron_cuts);
      LeptonFinder dressed_muons(bare_muons, photons, 0.1, muon_cuts);
      declare(dressed_electrons, "electrons");
      declare(dressed_muons, "muons");
      declare(photons, "photons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      book(_h_xsec_exp, 1, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve dressed leptons, sorted by pT
      Particles muons = apply<LeptonFinder>(event, "muons").particles();
      Particles electrons = apply<LeptonFinder>(event, "electrons").particles();

      Particles photons =
          apply<PromptFinalState>(event, "photons").particles(Cuts::abseta < 2.5 && Cuts::pT > 20 * GeV);
      idiscardIfAnyDeltaRLess(photons, muons, 0.5);
      idiscardIfAnyDeltaRLess(photons, electrons, 0.5);
      idiscardIfAnyDeltaRLess(electrons, muons, 0.5);

      // Retrieve Etmiss variable
      const double Etmiss = apply<MissingMomentum>(event, "MET").missingEt();

      if ((muons.size() < 1) || (electrons.size() < 1))  vetoEvent;
      if (photons.size() != 1)  vetoEvent;
      Particle muon = muons.at(0);
      Particle electron = electrons.at(0);
      if (muon.pid() * electron.pid() > 0)  vetoEvent;
      Particle photon = photons.at(0);
      FourMomentum LL = muon.momentum() + electron.momentum();
      if (LL.mass() < 10 * GeV)  vetoEvent;
      if (LL.pt() < 15 * GeV)    vetoEvent;

      FourMomentum EtMiss = apply<MissingMomentum>(event, "MET").missingMomentum();
      const double dphi = deltaPhi(LL, EtMiss);
      const double mT = sqrt(2 * LL.pT() * EtMiss.pT() * (1 - cos(dphi)));

      if (mT < 10*GeV)      vetoEvent;
      if (Etmiss < 20*GeV)  vetoEvent;
      _h_xsec_exp->fill("Results"s);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double norm = crossSection() / femtobarn / sumOfWeights();

      scale(_h_xsec_exp, norm);
    }

    BinnedHistoPtr<string> _h_xsec_exp;
  };

  RIVET_DECLARE_PLUGIN(CMS_2023_I2709669);

}  // namespace Rivet
