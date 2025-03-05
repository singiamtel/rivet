// -*- C++ -*-
#ifndef RIVET_DISLepton_HH
#define RIVET_DISLepton_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/UndressBeamLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {

  /// Enum for available DIS incoming-to-outgoing lepton mode
  enum class DISMode { L2L, L2NU, NU2L, NU2NU };

  /// @brief Get the incoming and outgoing leptons in a DIS event.
  class DISLepton : public FinalState {
  public:

    /// @name Constructors.
    /// @{

    /// Constructor with optional cuts first
    DISLepton(const Cut& cuts=Cuts::OPEN,
	      LeptonReco lreco=LeptonReco::ALL, ObjOrdering lsort=ObjOrdering::ENERGY,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0, DISMode dismode=DISMode::L2L)
      : _isolDR(isolDR), _lsort(lsort), _lreco(lreco), _dismode(dismode)
    {
      setName("DISLepton");
      declare(HadronicFinalState(), "IFS");

      // Beam undressing
      if (beamundresstheta > 0.0) {
        declare(UndressBeamLeptons(beamundresstheta), "Beam");
      } else {
        declare(Beam(), "Beam");
      }

      // Lepton reco mode
      switch (_lreco) {
      case LeptonReco::ALL:
        declare(FinalState(cuts), "LFS");
        break;
      case LeptonReco::ALL_DRESSED:
        declare(LeptonFinder(FinalState(), dressDR, cuts), "LFS");
        break;
      case LeptonReco::PROMPT_BARE:
        declare(PromptFinalState(cuts), "LFS");
        break;
      case LeptonReco::PROMPT_DRESSED:
        declare(LeptonFinder(PromptFinalState(), dressDR, cuts), "LFS");
        break;
      }

      // Identify the non-outgoing-lepton part of the event
      VetoedFinalState remainingFS;
      remainingFS.addVetoOnThisFinalState(*this);
      declare(remainingFS, "RFS");
    }


    /// Constructor without lepton-ordering spec, requiring cuts
    DISLepton(Cut& cuts, LeptonReco lreco=LeptonReco::ALL,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0, DISMode dismode=DISMode::L2L)
      : DISLepton(cuts, lreco, ObjOrdering::ENERGY, beamundresstheta, isolDR, dressDR, dismode)
    {  }

    /// Constructor without cuts, requiring lepton reco spec
    DISLepton(LeptonReco lreco, ObjOrdering lsort=ObjOrdering::ENERGY,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0, DISMode dismode=DISMode::L2L)
      : DISLepton(Cuts::OPEN, lreco, lsort, beamundresstheta, isolDR, dressDR, dismode)
    {  }

    /// Constructor without cuts or lepton-ordering spec, requiring lepton reco spec
    DISLepton(LeptonReco lreco,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0, DISMode dismode=DISMode::L2L)
      : DISLepton(Cuts::OPEN, lreco, ObjOrdering::ENERGY, beamundresstheta, isolDR, dressDR, dismode)
    {  }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(DISLepton);

    /// @}

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual CmpState compare(const Projection& p) const;


  public:

    /// The incoming lepton
    const Particle& in() const { return _incoming; }

    /// The outgoing lepton
    const Particle& out() const { return _outgoing; }

    /// Sign of the incoming lepton pz component
    int pzSign() const { return sign(_incoming.pz()); }

    /// Lepton reconstruction mode
    LeptonReco reconstructionMode() const { return _lreco; }


    /// Access to the particles other than outgoing leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;

    /// Clear the projection
    void clear() { _theParticles.clear(); }


  protected:

    /// The incoming lepton
    Particle _incoming;

    /// The outgoing lepton
    Particle _outgoing;

    /// An isolation cone around the lepton
    double _isolDR;

    /// How to sort leptons
    ObjOrdering _lsort;

    /// The lepton reconstruction mode
    LeptonReco _lreco;

    /// The incoming-to-outgoing lepton mode
    DISMode _dismode;

  };


}

#endif
