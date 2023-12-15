// -*- C++ -*-
#ifndef RIVET_USERCENTESTIMATE_HH
#define RIVET_USERCENTESTIMATE_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"

namespace Rivet {


  class UserCentEstimate: public SingleValueProjection {
  public:

    using SingleValueProjection::operator=;

    UserCentEstimate() {
      setName("UserCentEstimate");
      declare(HepMCHeavyIon(), "HepMC");
    }

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(UserCentEstimate);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    void project(const Event& e) {
      clear();
      setValue(apply<HepMCHeavyIon>(e, "HepMC").user_cent_estimate());
    }

    CmpState compare(const Projection&) const {
      return CmpState::EQ;
    }

  };


}

#endif
