#include "Rivet/Tools/JetUtils.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/FastJets.hh"

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

  Jet& itrimJetsFrac(Jet &jet, const double frac){
    // this keeps all (ghost-associated) tags
    // this approach might be incorrect
    // alternatively, could discard all tags, which might be incorrect as well
    // but: discarding tags is easily done by user, so let's keep everything
    const double ptcut = jet.pT()*frac;
    vector<Particle> preservedParts;
    FourMomentum fourmom;
    for (const Particle &jconstit : jet.constituents()){
      if (jconstit.pT() > ptcut){
        preservedParts.push_back(jconstit);
        fourmom += jconstit;
      }
    }
    jet.setState(fourmom, preservedParts, jet.tags());
    return jet;
  }

  PseudoJet& ifilterPseudoJets(PseudoJet &pj, const fastjet::Filter &filter){
    pj = filter(pj);
    return pj;
  }

}
