// -*- C++ -*-
#ifndef RIVET_DISKinematics_HH
#define RIVET_DISKinematics_HH

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// Type of DIS boost to apply
  enum class DISFrame { HCM, XCM, BREIT, LAB };

  
  /// @brief Get the DIS kinematic variables and relevant boosts for an event.
  class DISKinematics : public Projection {
  public:

    /// The default constructor.
    DISKinematics(const DISLepton& lepton=DISLepton())
      : _theQ2(-1.0), _theW2(-1.0), _theX(-1.0), _theY(-1.0), _theS(-1.0), _theGH(-1.0)
    {
      setName("DISKinematics");
      //addPdgIdPair(ANY, hadid);
      declare(Beam(), "Beam");
      declare(lepton, "Lepton");
    }

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(DISKinematics);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual CmpState compare(const Projection& p) const;


  public:

    /// @brief The \f$Q^2\f$, i.e. the virtuality of the DIS photon.
    ///
    /// Calculated by identifying the virtual photon momentum as the
    /// difference between incoming and outgoing scattered-lepton momenta.
    double Q2() const { return _theQ2; }

    /// The \f$W^2\f$, i.e. the squared CoM energy of the DIS-photon+hadron system.
    double W2() const { return _theW2; }

    /// The Bjorken \f$x\f$.
    double x() const { return _theX; }

    /// The inelasticity \f$y\f$
    double y() const { return _theY; }

    /// The centre of mass energy \f$s\f$
    double s() const { return _theS; }

    /// The angle \$f\gamma_{had}\f$
    double gammahad() const { return _theGH; }



    /// The LorentzRotation needed to boost a particle to the hadronic CM frame.
    const LorentzTransform& boostHCM() const {
      return _hcm;
    }

    /// The LorentzRotation needed to boost a particle to the hadronic Breit frame.
    const LorentzTransform& boostBreit() const {
      return _breit;
    }

    /// The incoming hadron beam particle
    const Particle& beamHadron() const {
      return _inHadron;
    }

    /// The incoming lepton beam particle
    const Particle& beamLepton() const {
      return _inLepton;
    }

    /// The scattered DIS lepton
    const Particle& scatteredLepton() const {
      return _outLepton;
    }

    /// @brief 1/-1 multiplier indicating (respectively) whether the event has conventional orientation or not
    ///
    /// Conventional DIS orientation has the hadron travelling in the +z direction
    int orientation() const {
      return sign(_inHadron.pz());
    }


  protected:

    /// The \f$Q^2\f$.
    double _theQ2;

    /// The \f$W^2\f$.
    double _theW2;

    /// The Bjorken \f$x\f$.
    double _theX;

    /// The Inelasticity \f$y\f$
    double _theY;

    /// The centre of mass energy \f$s\f$
    double _theS;
   /// The angle \f$\gamma_{had}\f$
    double _theGH;


    /// Incoming and outgoing DIS particles
    Particle _inHadron, _inLepton, _outLepton;

    /// The LorentzRotation needed to boost a particle to the hadronic CM frame.
    LorentzTransform _hcm;

    /// The LorentzRotation needed to boost a particle to the hadronic Breit frame.
    LorentzTransform _breit;

  };


}

#endif
