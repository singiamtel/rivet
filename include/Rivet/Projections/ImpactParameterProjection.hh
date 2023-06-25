// -*- C++ -*-
#ifndef RIVET_IMPACTPARAMETERPROJECTION_HH
#define RIVET_IMPACTPARAMETERPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"


namespace Rivet {

  class ImpactParameterProjection: public SingleValueProjection {
  public:

    ImpactParameterProjection() {
      setName("ImpactParameterProjection");
      declare(HepMCHeavyIon(), "HepMC");
    }

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(ImpactParameterProjection);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    void project(const Event& e) {
      clear();
      set(apply<HepMCHeavyIon>(e, "HepMC").impact_parameter());
    }

    CmpState compare(const Projection& p) const {
      return CmpState::EQ;
    }

  };

}

#endif
