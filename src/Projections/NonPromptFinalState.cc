// -*- C++ -*-
#include "Rivet/Projections/NonPromptFinalState.hh"

namespace Rivet {


  NonPromptFinalState::NonPromptFinalState(const FinalState& fsp, TauDecaysAs taudecays, MuDecaysAs mudecays)
    : _mudecays(mudecays == MuDecaysAs::NONPROMPT),
      _taudecays(taudecays == TauDecaysAs::NONPROMPT)
  {
    setName("NonPromptFinalState");
    declare(fsp, "FS");
  }


  NonPromptFinalState::NonPromptFinalState(const Cut& c, TauDecaysAs taudecays, MuDecaysAs mudecays)
    : _mudecays(mudecays == MuDecaysAs::NONPROMPT),
      _taudecays(taudecays == TauDecaysAs::NONPROMPT)
  {
    setName("NonPromptFinalState");
    declare(FinalState(c), "FS");
  }


  CmpState NonPromptFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != CmpState::EQ) return fscmp;
    const NonPromptFinalState& other = dynamic_cast<const NonPromptFinalState&>(p);
    return cmp(_mudecays, other._mudecays) || cmp(_taudecays, other._taudecays);
  }


  void NonPromptFinalState::project(const Event& e) {
    _theParticles.clear();

    const Particles& particles = apply<FinalState>(e, "FS").particles();
    for (const Particle& p : particles)
      if (!isPrompt(p, !_taudecays, !_mudecays)) _theParticles.push_back(p);
    MSG_DEBUG("Number of final state particles from hadron decays = " << _theParticles.size());

    if (getLog().isActive(Log::TRACE)) {
      for (const Particle& p : _theParticles)
        MSG_TRACE("Selected: " << p.pid() << ", charge = " << p.charge());
    }
  }


}
