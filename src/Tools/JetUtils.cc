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


}
