// -*- C++ -*-
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/NonPromptFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  // Separate-FS constructor
  /// @note The first two args were swapped in version 3.2.0!
  LeptonFinder::LeptonFinder(const FinalState& leptonfs, const FinalState& photonfs,
                             double dRdress, const Cut& cut, DressingType dressing)
    : FinalState(cut),
      _dRdress(dRdress),
      _dressMode(dressing)
  {
    setName("LeptonFinder");

    IdentifiedFinalState leptonfs2(leptonfs); // make sure it really is only leptons
    leptonfs2.acceptIdPairs({PID::ELECTRON, PID::MUON, PID::TAU}); //< usually no final-state taus...
    IdentifiedFinalState photonfs2(photonfs, PID::PHOTON); // make sure it really is only photons
    declare(leptonfs2, "Leptons");
    declare(photonfs2, "Photons");

    // Set up FJ clustering
    if (_dressMode == DressingType::AKT) {
      MergedFinalState mergedfs(photonfs, leptonfs);
      FastJets leptonjets(mergedfs, JetAlg::ANTIKT, dRdress);
      declare(leptonjets, "LeptonJets");
    }
  }


  // No-FS version
  LeptonFinder::LeptonFinder(double dRdress, const Cut& cut,
                             LeptonOrigin whichleptons, PhotonOrigin whichphotons,
                             TauDecaysAs tauDecays, MuDecaysAs muDecays,
                             DressingType dressing)
    : FinalState(cut),
      _dRdress(dRdress),
      _dressMode(dressing)
  {
    setName("LeptonFinder");

    // Find leptons -- specialising to prompt if requested
    IdentifiedFinalState leptonfs;
    leptonfs.acceptIdPairs({PID::ELECTRON, PID::MUON, PID::TAU}); //< usually no final-state taus...
    if (whichleptons == LeptonOrigin::NODECAY) {
      declare(PromptFinalState(leptonfs, tauDecays, muDecays), "Leptons");
    } else if (whichleptons == LeptonOrigin::ALL) {
      declare(leptonfs, "Leptons");
    } else if (whichleptons == LeptonOrigin::DECAY) {
      declare(NonPromptFinalState(leptonfs, tauDecays, muDecays), "Leptons");
    } else
      throw LogicError("LeptonFinder cannot retrieve leptons when no origins are accepted");

    // Find photons -- specialising to prompt if requested
    /// @todo Generalise to allow other clustering particles, e.g. e+e- (don't double-count)
    IdentifiedFinalState photonfs(PID::PHOTON);
    if (whichphotons == PhotonOrigin::NODECAY) {
      declare(PromptFinalState(photonfs, tauDecays, muDecays), "Photons");
    } else if (whichphotons == PhotonOrigin::NONE) {
      declare(FinalState(photonfs, Cuts::abspid != PID::PHOTON), "Photons");
    } else {
      declare(photonfs, "Photons");
    }

    MergedFinalState mergedfs(photonfs, leptonfs);
    declare(mergedfs, "Constituents");

    // Set up FJ clustering
    if (_dressMode == DressingType::AKT) {
      FastJets leptonjets(mergedfs, JetAlg::ANTIKT, dRdress);
      declare(leptonjets, "LeptonJets");
    }
  }


  CmpState LeptonFinder::compare(const Projection& p) const {
    // Do the fast comparison of local stuff first
    const LeptonFinder& other = dynamic_cast<const LeptonFinder&>(p);
    CmpState localcmp = cmp(_dRdress, other._dRdress) || cmp(_dressMode, other._dressMode);
    if (localcmp != CmpState::EQ) return localcmp;

    // Compare the two projs as final states (for the cuts)
    CmpState fscmp = FinalState::compare(other);
    if (fscmp != CmpState::EQ) return fscmp;

    const PCmp phcmp = mkNamedPCmp(p, "Photons");
    if (phcmp != CmpState::EQ) return phcmp;

    const PCmp sigcmp = mkNamedPCmp(p, "Leptons");
    if (sigcmp != CmpState::EQ) return sigcmp;

    /// @todo Add a mode for the flavoured dressing
    if (_dressMode == DressingType::AKT) {
      const PCmp ljcmp = mkNamedPCmp(p, "LeptonJets");
      if (ljcmp != CmpState::EQ) return ljcmp;
    }

    return CmpState::NEQ;
  }


  void LeptonFinder::project(const Event& e) {
    _theParticles.clear();

    // Get bare leptons
    const Particles& bareleptons = apply<FinalState>(e, "Leptons").particles();
    MSG_DEBUG("Number of bare leptons: " << bareleptons.size());
    if (bareleptons.empty()) return;

    // Initialise DL collection with bare leptons
    Particles allClusteredLeptons;
    allClusteredLeptons.reserve(bareleptons.size());

    if (_dRdress <= 0 || _dressMode == DressingType::CONE) {

      // If the radius is 0 or negative, don't even attempt to cluster
      for (const Particle& bl : bareleptons) {
        Particle dl(bl.pid(), bl.momentum(), bl.genParticle(), bl.origin());
        dl.setConstituents({bl});
        allClusteredLeptons += dl;
      }

    }

    // Decide how to perform the clustering
    switch (_dressMode) {

      case DressingType::CONE:

        {
          // Match each photon to its closest charged lepton within the dR cone
          const FinalState& photons = apply<FinalState>(e, "Photons");
          for (const Particle& photon : photons.particles()) {
            double dRmin = _dRdress;
            int idx = -1;
            for (size_t i = 0; i < bareleptons.size(); ++i) {
              const Particle& bl = bareleptons[i];
              // Only cluster photons around *charged* signal particles
              if (bl.charge3() == 0) continue;
              // Find the closest lepton
              const double dR = deltaR(bl, photon);
              if (dR < dRmin) {
                dRmin = dR;
                idx = i;
              }
            }
            // Escape if no lepton found within the dRdress range
            if (idx < 0) continue;

            // Attach the photon to the closest in-range lepton
            Particle& dl = allClusteredLeptons[idx];
            MSG_DEBUG("Adding photon " << photon << " to dressed lepton #" << idx << ": " << dl);
            dl.addConstituent(photon, true);
          }
        }
        break;


      // case DressingType::AKT:
      case DressingType::CLUSTER:

        {
          if (allClusteredLeptons.size()) {
            allClusteredLeptons.clear();
            allClusteredLeptons.reserve(bareleptons.size());
          }
          const Jets& lepjets = apply<JetFinder>(e, "LeptonJets").jets();
          for (const Jet& lepjet : lepjets) {
            const Particles leps = sortByPt(lepjet.particles(isChargedLepton));
            if (leps.empty()) continue;
            Particles constituents = {leps[0]}; //< note no dressing for subleading leptons
            Particle dl(leps[0].pid(), leps[0].momentum(), leps[0].genParticle(), leps[0].origin());
            constituents += lepjet.particles(isPhoton);
            dl.setConstituents(constituents);
            allClusteredLeptons += dl;
          }
        }
        break;

      // case DressingType::FLAVKT:
      // case DressingType::FLAVCLUSTER:

        // {
        //   const Jets& lepjets = apply<JetFinder>(e, "LeptonJets").jets();
        //   for (const Jet& lepjet : lepjets) {
        //     const Particles leps = sortByPt(lepjet.particles(isChargedLepton));
        //     if (leps.empty()) continue;
        //     Particles constituents = {leps[0]}; //< note no dressing for subleading leptons
        //     Particle dl(leps[0].pid(), leps[0].momentum(), leps[0].genParticle(), leps[0].origin());
        //     constituents += lepjet.particles(isPhoton);
        //     dl.setConstituents(constituents);
        //     allClusteredLeptons += dl;
        //   }
        // }
        // break;

    }


    // Fill the canonical particles collection with the composite DL Particles
    for (const Particle& lepton : allClusteredLeptons) {
      const bool acc = accept(lepton);
      MSG_TRACE("Clustered lepton " << lepton
                << " with constituents = " << lepton.constituents()
                << ", cut-pass = " << std::boolalpha << acc);
      if (acc) _theParticles.push_back(lepton);
    }
    MSG_DEBUG("#dressed leptons = " << allClusteredLeptons.size()
              << " -> " << _theParticles.size() << " after cuts");

  }


}
