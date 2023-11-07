// -*- C++ -*-
#ifndef RIVET_InvisibleFinalState_HH
#define RIVET_InvisibleFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  enum class OnlyPrompt { YES, NO };
  
  
  /// @brief Final state modifier excluding particles which are experimentally visible
  class InvisibleFinalState : public FinalState {
  public:

    /// @name Constructors
    /// @{

    /// Constructor with specific FinalState.
    InvisibleFinalState(OnlyPrompt requirepromptness=OnlyPrompt::NO,
			TauDecaysAs taudecays=TauDecaysAs::NONPROMPT,
			MuDecaysAs mudecays=MuDecaysAs::NONPROMPT)
      : _requirePromptness(requirepromptness == OnlyPrompt::YES),
        _taudecays(taudecays == TauDecaysAs::PROMPT),
        _mudecays(mudecays == MuDecaysAs::PROMPT)
    {
      setName("InvisibleFinalState");
      declare(FinalState(), "FS");
    }

    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(InvisibleFinalState);

    /// @}

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Require accepted particles to be prompt
    void requirePromptness(bool acc=true,
			   TauDecaysAs taudecays=TauDecaysAs::NONPROMPT,
			   MuDecaysAs mudecays=MuDecaysAs::NONPROMPT) {
      _requirePromptness = acc;
      _taudecays = (taudecays == TauDecaysAs::PROMPT);
      _mudecays = (mudecays == MuDecaysAs::PROMPT);
    }

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  protected:

    bool _requirePromptness, _taudecays, _mudecays;
  };


}

#endif
