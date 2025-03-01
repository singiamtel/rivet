// -*- C++ -*-
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  void UnstableParticles::project(const Event& e) {
    _theParticles.clear();

    /// @todo Replace PID veto list with PID:: functions?
    vector<PdgId> vetoIds;
    vetoIds += PID::PHOTON; // status 2 photons don't count!
    vetoIds += 110; vetoIds += 990; vetoIds += 9990; // Reggeons

    for (ConstGenParticlePtr p : HepMCUtils::particles(e.genEvent())) {

      // First select on GenParticle info/object only
      bool passed =
        (p->status() == 1 || //< Stable particle... huh?!?
         (p->status() == 2 && !contains(vetoIds, abs(p->pdg_id())))) && //< Unvetoed unstable
        !PID::isParton(p->pdg_id()) && //< Always veto partons, regardless of status
        p->status() != 4; //< Filter beam particles (necessary?)

      // Avoid double-counting by only using the last item in a decay chain
      /// @todo Replace with isLastWith(pid)?
      /// @todo Generalise to allow selection of first or last particles in duplicate chains
      if (passed) {
        if (p->end_vertex()) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(p->end_vertex(), Relatives::CHILDREN)) {
            if (p->pdg_id() == pp->pdg_id() && pp->status() == 2) {
              passed = false;
              break;
            }
          }
        }
      }

      // Apply user cuts on the Particle wrapper and add to result if successful
      if (passed) {
        const Particle rp(p);
        if (_cuts->accept(rp)) _theParticles.push_back(std::move(rp));
      }

      // Log parents and children (hence not short-circuiting vetos above)
      if (getLog().isActive(Log::TRACE)) {
        MSG_TRACE("ID = " << p->pdg_id()
                  << ", status = " << p->status()
                  << ", pT = " << p->momentum().perp()
                  << ", eta = " << p->momentum().eta()
                  << ": result = " << std::boolalpha << passed);
        if (p->production_vertex()) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(p->production_vertex(), Relatives::PARENTS))
            MSG_TRACE("  parent ID = " << pp->pdg_id());
        }
        if (p->end_vertex()) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(p->end_vertex(), Relatives::CHILDREN))
            MSG_TRACE("  child ID  = " << pp->pdg_id());
        }
      }
    }
    MSG_DEBUG("Number of unstable final-state particles = " << _theParticles.size());
  }


}
