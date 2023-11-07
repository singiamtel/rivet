// -*- C++ -*-
#ifndef RIVET_NonPromptFinalState_HH
#define RIVET_NonPromptFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.fhh"

namespace Rivet {

  
  /// @brief Find final state particles NOT directly connected to the hard process.
  ///
  /// See PromptFinalState for details.
  class NonPromptFinalState : public FinalState {
  public:

    /// @name Constructors
    /// @{

    // Constructor from a final state.
    NonPromptFinalState(const FinalState& fsp,
			TauDecaysAs taudecays=TauDecaysAs::PROMPT,
			MuDecaysAs mudecays=MuDecaysAs::PROMPT);

    /// Constructor from a Cut (and implicit general FS).
    NonPromptFinalState(const Cut& c,
			TauDecaysAs taudecays=TauDecaysAs::PROMPT,
			MuDecaysAs mudecays=MuDecaysAs::PROMPT);

    // /// Constructor from a Cut and optional FinalState.
    // NonPromptFinalState(const Cut& c, const FinalState& fsp=FinalState(),
    // 			TauDecaysAs taudecays=TauDecaysAs::PROMPT,
    // 			MuDecaysAs mudecays=MuDecaysAs::PROMPT);

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(NonPromptFinalState);

    /// @}

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Treat particles from decays of prompt muons as non-prompt?
    void acceptMuonDecays(bool acc=true) { _mudecays = acc; }
    /// Treat particles from decays of prompt taus as non-prompt?
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
