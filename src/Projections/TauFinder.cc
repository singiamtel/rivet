// -*- C++ -*-
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  void TauFinder::project(const Event& e) {
    _theParticles.clear();

    const Particles& ufs = apply<UnstableParticles>(e, "UFS").particles();
    for (const Particle& p : ufs) {
      // Taus only, obv
      if (p.abspid() != PID::TAU) continue;

      // Discard if wrong decay mode
      if ((_decmode == TauDecay::LEPTONIC && isHadronic(p)) ||
	  (_decmode == TauDecay::HADRONIC && isLeptonic(p)))
        continue;

      // Discard if wrong origin
      if ((_origin == LeptonOrigin::NODECAY && !p.isPrompt()) ||
          (_origin == LeptonOrigin::DECAY && p.isPrompt()))
        continue;
      if (_origin == LeptonOrigin::NONE)
        throw LogicError("TauFinder cannot retrieve taus when no origins are accepted");

      // If not discarded
      _theParticles.push_back(p);
    }
  }


  CmpState TauFinder::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "UFS");
    if (fscmp != CmpState::EQ) return fscmp;

    const TauFinder& other = dynamic_cast<const TauFinder&>(p);
    return cmp(_decmode, other._decmode) || cmp(_origin, other._origin);
  }


}
