// -*- C++ -*-
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  DileptonFinder::DileptonFinder(double masstarget, double dRmax,
                                 const Cut& lcuts, const Cut& llcuts,
                                 LeptonOrigin whichleptons, PhotonOrigin whichphotons,
                                 TauDecaysAs tauDecays, MuDecaysAs muDecays,
                                 DressingType dressing)
    : FinalState(llcuts), _masstarget(masstarget)
  {
    setName("DileptonFinder");

    LeptonFinder leptons(lcuts, dRmax, whichleptons, whichphotons, tauDecays, muDecays, dressing);
    declare(leptons, "Leptons");

    // Identify the non-Z part of the event
    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    declare(remainingFS, "RFS");
  }


  DileptonFinder::DileptonFinder(const FinalState& allfs,
                                 double masstarget, double dRmax,
                                 const Cut& lcuts, const Cut& llcuts,
                                 DressingType dressing)
    : FinalState(llcuts), _masstarget(masstarget)
  {
    setName("DileptonFinder");

    LeptonFinder leptons(allfs, lcuts, dRmax, dressing);
    declare(leptons, "Leptons");

    // Identify the non-Z part of the event
    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    declare(remainingFS, "RFS");
  }


  DileptonFinder::DileptonFinder(const FinalState& leptonfs, const FinalState& photonfs,
                                 double masstarget, double dRmax,
                                 const Cut& lcuts, const Cut& llcuts,
                                 DressingType dressing)
    : FinalState(llcuts), _masstarget(masstarget)
  {
    setName("DileptonFinder");

    LeptonFinder leptons(photonfs, leptonfs, ///< @todo Swap!!
                         dRmax, lcuts, dressing);
    declare(leptons, "Leptons");

    // Identify the non-Z part of the event
    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    declare(remainingFS, "RFS");
  }


  /////////////////////////////////////////////////////


  const Particles& DileptonFinder::constituents() const {
    static const Particles NO_PARTICLES;
    return empty() ? NO_PARTICLES : boson().constituents();
  }


  const VetoedFinalState& DileptonFinder::remainingFinalState() const {
    return getProjection<VetoedFinalState>("RFS");
  }


  CmpState DileptonFinder::compare(const Projection& p) const {
    const DileptonFinder& other = dynamic_cast<const DileptonFinder&>(p);
    CmpState mcmp = cmp(_masstarget, other._masstarget);
    if (mcmp != CmpState::EQ) return mcmp;

    // Compare cuts via FS
    CmpState fscmp = FinalState::compare(other);
    if (fscmp != CmpState::EQ) return fscmp;

    // Compare lepton finders
    PCmp lepcmp = mkNamedPCmp(p, "Leptons");
    if (lepcmp != CmpState::EQ) return lepcmp;

    return CmpState::NEQ;
  }


  void DileptonFinder::project(const Event& e) {
    clear();

    // Get leptons and find an acceptable invariant mass OSSF pair
    const Particles& leptons = apply<LeptonFinder>(e, "Leptons").particles();
    const Particles leptonsp = select(leptons, Cuts::charge > 0);
    const Particles leptonsm = select(leptons, Cuts::charge < 0);
    pair<int,int> ij = closestMatchIndices(leptonsp, leptonsm, Kin::mass, 91.2*GeV);
    if (ij.first < 0 || ij.second < 0) {
      MSG_TRACE("No acceptable inv-mass l+l- pairs found");
      return;
    }

    // InvMassFinalState imfs({_pid, -_pid}, _minmass, _maxmass, _masstarget);
    // imfs.calc(leptons.particles());
    // if (imfs.particlePairs().empty()) {
    //   MSG_TRACE("No acceptable inv-mass lepton/antilepton pairs found");
    //   return;
    // }

    // Assemble a pseudo-Z particle
    const Particle& pp(leptonsp[ij.first]), pm(leptonsm[ij.second]);
    const FourMomentum pZ = pp.momentum() + pm.momentum();
    assert(pp.charge3() + pm.charge3() == 0);
    Particle z(PID::Z0BOSON, pZ);
    MSG_DEBUG(z << " reconstructed from: " << pp << " + " << pm);

    // Add (dressed) lepton constituents to the Z
    // Keep the dressed leptons found by the DileptonFinder
    MSG_TRACE("l+ = " << pp.constituents());
    MSG_TRACE("l- = " << pm.constituents());
    z.addConstituent(pp);
    z.addConstituent(pm);
    MSG_DEBUG("Number of stored raw Z constituents = "
              << z.rawConstituents().size() << "  " << z.rawConstituents());

    // Register the completed Z
    _theParticles.push_back(z);

    // Apply the ll filtering cuts (registered at ParticleFinder level)
    iselect(_theParticles, _cuts);
  }


}
