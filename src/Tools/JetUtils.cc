#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  Jets& iselect(Jets& jets, const Cut& c) {
    if (c == Cuts::OPEN) return jets;
    // return iselect(jets, *c);
    return iselect(jets, [&](const Jet& j){return c->accept(j);});
  }


  Jets& idiscard(Jets& jets, const Cut& c) {
    if (c == Cuts::OPEN) { jets.clear(); return jets; }
    // return idiscard(jets, *c);
    return idiscard(jets, [&](const Jet& j){return c->accept(j);});
  }

  Jets trimJetsFrac(const PseudoJets& jetsIn, const double frac){
    Jets jetsOut;
    for (const PseudoJet & pj : jetsIn){
      const double ptcut = pj.pt()*frac;
      // TODO: Ugly -> I'd like to use Jets instead of Particles
      // Will need changes to jet class.
      vector<Particle> preservedJets;
      FourMomentum fourmom;
      for (const PseudoJet & pjconstit : pj.constituents()){
        if (pjconstit.pt() > ptcut){
          preservedJets.push_back(Particle(0, Jet(pjconstit).mom()));
          fourmom += Jet(pjconstit);
        }
      }
      jetsOut.emplace_back(fourmom, preservedJets);
    }
    return jetsOut;
  }

}
