#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// @brief dileptonic e-mu ttbar
  class ATLAS_2016_I1468168 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1468168);

    void init() {
      // Eta ranges
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT >= 1.0*MeV;
      // Lepton cuts
      Cut lep_cuts = Cuts::abseta < 2.5 && Cuts::pT >= 25.0*GeV;

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      LeptonFinder dressedelectrons(electrons, photons, 0.1, lep_cuts);
      declare(dressedelectrons, "DressedElectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      LeptonFinder dressedmuons(muons, photons, 0.1, lep_cuts);
      declare(dressedmuons, "DressedMuons");

      book(_h , 2, 1 ,1);
    }


    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      const size_t num_es = apply<LeptonFinder>(event, "DressedElectrons").dressedLeptons().size();
      const size_t num_mus = apply<LeptonFinder>(event, "DressedMuons").dressedLeptons().size();

      // Evaluate basic event selection
      const bool pass_emu = num_es == 1 && num_mus == 1;
      if (!pass_emu) vetoEvent;

      // Fill histogram to measure the event acceptance
      _h->fill(13000);
    }


    void finalize() {
      // Normalize to cross-section
      scale(_h, crossSection()/picobarn / sumOfWeights());
    }


   private:

    BinnedHistoPtr<int> _h;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1468168);

}
