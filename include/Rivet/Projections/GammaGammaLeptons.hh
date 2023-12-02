// -*- C++ -*-
#ifndef RIVET_GammaGammaLeptons_HH
#define RIVET_GammaGammaLeptons_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/UndressBeamLeptons.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Get the incoming and outgoing leptons in a gamma gamma collision event in e+e-
  class GammaGammaLeptons : public Projection {
  public:

    /// @name Constructors.
    /// @{

    /// Constructor with optional cuts first
    GammaGammaLeptons(const Cut& cuts=Cuts::OPEN,
		      LeptonReco lreco=LeptonReco::ALL, ObjOrdering lsort=ObjOrdering::ENERGY,
		      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0)
      : _isolDR(isolDR), _lsort(lsort)
    {
      setName("GammaGammaLeptons");
      declare(HadronicFinalState(), "IFS");

      // Beam undressing
      if (beamundresstheta > 0.0) {
        declare(UndressBeamLeptons(beamundresstheta), "Beam");
      } else {
        declare(Beam(), "Beam");
      }

      // Lepton reco mode
      switch (lreco) {
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
    }


    /// Constructor without lepton-ordering spec, requiring cuts
    GammaGammaLeptons(Cut& cuts, LeptonReco lreco=LeptonReco::ALL,
		      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0)
      : GammaGammaLeptons(cuts, lreco, ObjOrdering::ENERGY, beamundresstheta, isolDR, dressDR)
    {  }

    /// Constructor without cuts, requiring lepton reco spec
    GammaGammaLeptons(LeptonReco lreco, ObjOrdering lsort=ObjOrdering::ENERGY,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0)
      : GammaGammaLeptons(Cuts::OPEN, lreco, lsort, beamundresstheta, isolDR, dressDR)
    {  }

    /// Constructor without cuts or lepton-ordering spec, requiring lepton reco spec
    GammaGammaLeptons(LeptonReco lreco,
	      double beamundresstheta=0.0, double isolDR=0.0, double dressDR=0.0)
      : GammaGammaLeptons(Cuts::OPEN, lreco, ObjOrdering::ENERGY, beamundresstheta, isolDR, dressDR)
    {  }


    /// Constructor from other constructors
    GammaGammaLeptons(const FinalState& leptoncandidates,
		      const Beam& beamproj=Beam(),
		      const FinalState& isolationfs=FinalState(),
		      double isolDR = 0.0, ObjOrdering lsort=ObjOrdering::ENERGY)
      : _isolDR(isolDR), _lsort(lsort)
    {
      declare(leptoncandidates, "LFS");
      declare(isolationfs, "IFS");
      declare(beamproj, "Beam");
    }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(GammaGammaLeptons);

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
    const ParticlePair & in() const { return _incoming; }

    /// The outgoing lepton
    const ParticlePair & out() const { return _outgoing; }


  protected:

    /// The incoming leptons
    ParticlePair _incoming;

    /// The outgoing leptons
    ParticlePair _outgoing;

    /// An isolation cut around the lepton
    double _isolDR;

    /// How to sort leptons
    ObjOrdering _lsort;

  };


}

#endif
