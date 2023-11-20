// -*- C++ -*-
#ifndef RIVET_DISFinalState_HH
#define RIVET_DISFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Final state particles boosted to the hadronic center of mass system.
  ///
  /// NB. The DIS scattered lepton is not included in the final state particles.
  class DISFinalState: public FinalState {
  public:

    /// @name Constructors
    /// @{

    /// @brief Constructor with explicit FinalState
    DISFinalState(const FinalState& fs, DISFrame boosttype)
      : _boosttype(boosttype)
    {
      setName("DISFinalState");
      declare(fs, "FS");
      declare(DISKinematics(), "Kinematics");
    }


    /// @brief Constructor with optional FinalState
    DISFinalState(DISFrame boosttype, const FinalState& fs=FinalState())
      : DISFinalState(fs, boosttype)
    {    }


    /// @brief Constructor with explicit cuts to define final-state particles
    ///
    /// @note The cuts will be applied *before* the boost, e.g. to express detector acceptance.
    ///
    /// @todo Add a second optional Cut argument for post-boost cuts.
    DISFinalState(const Cut& c, DISFrame boosttype)
      : DISFinalState(FinalState(c), boosttype)
    {    }


    /// @brief Constructor with explicit cuts to define final-state particles
    ///
    /// @note The cuts will be applied *before* the boost, e.g. to express detector acceptance.
    ///
    /// @todo Add a second optional Cut argument for post-boost cuts.
    ///
    /// @deprecated The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(DISFrame boosttype, const Cut& c)
      : DISFinalState(FinalState(c), boosttype)
    {    }


    // /// @brief Constructor with default FinalState
    // DISFinalState(BoostFrame boosttype)
    //   : DISFinalState(FinalState(), boosttype)
    // {    }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(DISFinalState);

    /// @}


    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Get the associated DISKinematics (to avoid needing a separate projection)
    const DISKinematics& kinematics() {
      return getProjection<DISKinematics>("Kinematics");
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    
    /// Compare projections.
    CmpState compare(const Projection& p) const {
      const DISFinalState& other = dynamic_cast<const DISFinalState&>(p);
      return mkNamedPCmp(p, "Kinematics") || mkNamedPCmp(p, "FS") || cmp(_boosttype, other._boosttype);
    }


  protected:

    DISFrame _boosttype;

  };


}

#endif
