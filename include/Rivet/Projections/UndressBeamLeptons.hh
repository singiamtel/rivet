// -*- C++ -*-
#ifndef RIVET_UndressBeamLeptons_HH
#define RIVET_UndressBeamLeptons_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Incoming lepton beams with collinear photons subtracted
  class UndressBeamLeptons : public Beam {
  public:

    using Beam::operator=;

    /// @brief Constructor with an angle as argument
    ///
    /// The momentum of any photon within this angle w.r.t. a charged
    /// lepton beam will be subtracted from the beam-lepton momentum.
    UndressBeamLeptons(double theta = 0.0): _thetamax(theta) {
      setName("UndressBeamLeptons");
      declare(FinalState(), "FS");
    }

    /// Clone on the heap
    RIVET_DEFAULT_PROJ_CLONE(UndressBeamLeptons);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Project on to the Event
    virtual void project(const Event& e);


  protected:

    /// Compare with other projections.
    virtual CmpState compare(const Projection & p) const;

    /// The beam particles in the current collision
    double _thetamax;

  };


}

#endif
