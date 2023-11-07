#ifndef RIVET_TauFinder_HH
#define RIVET_TauFinder_HH

#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// Enumerate types of tau decay
  enum class TauDecay {
    ANY = 0,
    ALL = 0,
    LEPTONIC,
    HADRONIC
  };


  
  /// @brief Convenience finder of unstable taus
  ///
  /// @todo Convert to a general ParticleFinder, but this has many knock-on
  /// effects and requires e.g. FastJets to also be generalised.
  class TauFinder : public FinalState {
  public:

    static bool isHadronic(const Particle& tau) {
      assert(tau.abspid() == PID::TAU);
      return any(tau.stableDescendants(), isHadron);
    }

    static bool isLeptonic(const Particle& tau) {
      return !isHadronic(tau);
    }


    TauFinder(TauDecay decaymode=TauDecay::ANY, const Cut& cut=Cuts::open()) {
      /// @todo What about directness/promptness?
      setName("TauFinder");
      _decmode = decaymode;
      declare(UnstableParticles(cut), "UFS");
    }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(TauFinder);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    const Particles& taus() const { return _theParticles; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare with other projections.
    virtual CmpState compare(const Projection& p) const;


  protected:

    /// The decay-mode enum
    TauDecay _decmode;

  };


  /// @todo Make this the canonical name in future?
  using Taus = TauFinder;


}


#endif
