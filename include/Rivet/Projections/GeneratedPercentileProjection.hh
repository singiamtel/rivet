// -*- C++ -*-
#ifndef RIVET_GENERATEDPERCENTILEPROJECTION_HH
#define RIVET_GENERATEDPERCENTILEPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {

  class GeneratedPercentileProjection: public SingleValueProjection {
  public:

  using SingleValueProjection::operator=;

  GeneratedPercentileProjection() {
    setName("GeneratedPercentileProjection");
    declare(HepMCHeavyIon(), "HepMC");
  }

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(GeneratedPercentileProjection);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    void project(const Event& e) {
      clear();
      setValue(apply<HepMCHeavyIon>(e, "HepMC").centrality());
    }

    CmpState compare(const Projection&) const {
      return CmpState::EQ;
    }

  };

}

#endif
