// -*- C++ -*-
#ifndef RIVET_PromptFinalState_HH
#define RIVET_PromptFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.fhh"

namespace Rivet {
  

  /// @brief Find final state particles directly connected to the hard process.
  ///
  /// The definition of "prompt" used in Rivet is that from high-scale physics, i.e.
  /// particles directly connected to the hard process in an interaction, regardless
  /// of realistic reconstructibility of displaced vertices, etc. By construction
  /// hadrons cannot be considered prompt as they will be colour connected to other
  /// parts of the event through non-perturbative effects: this projection can
  /// return electrons, muons, photons, and exotic particles which do not have a
  /// hadron in their post-hadronization ancestor chain. Flags exist to choose
  /// whether intermediate tau or muon decays invalidate a particle's promptness.
  ///
  /// @todo Decide how to treat brem photons off prompt leptons -- are they also prompt? "Decay" does not change the lepton PID...
  ///
  /// @todo Use enums for tau, mu, brem
  class PromptFinalState : public FinalState {
  public:

    /// @name Constructors
    /// @{

    /// Constructor without cuts
    PromptFinalState(TauDecaysAs taudecays=TauDecaysAs::NONPROMPT, MuDecaysAs mudecays=MuDecaysAs::NONPROMPT);

    /// Constructor from a Cut
    PromptFinalState(const Cut& c, TauDecaysAs taudecays=TauDecaysAs::NONPROMPT, MuDecaysAs mudecays=MuDecaysAs::NONPROMPT);

    // Constructor from a FinalState
    PromptFinalState(const FinalState& fsp, TauDecaysAs taudecays=TauDecaysAs::NONPROMPT, MuDecaysAs mudecays=MuDecaysAs::NONPROMPT);

    // /// Constructor from a Cut and optional FinalState.
    // PromptFinalState(const Cut& c, const FinalState& fsp=FinalState(), TauDecaysAs taudecays, MuDecaysAs mudecays);

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(PromptFinalState);

    /// @}

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Accept leptons from decays of prompt muons as themselves being prompt?
    void acceptMuonDecays(bool acc=true) { _mudecays = acc; }
    /// Accept leptons from decays of prompt taus as themselves being prompt?
    void acceptTauDecays(bool acc=true) { _taudecays = acc; }


    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  protected:

    bool _mudecays, _taudecays;

  };


}

#endif
