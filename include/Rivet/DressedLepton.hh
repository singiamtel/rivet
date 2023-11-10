// -*- C++ -*-
#ifndef RIVET_DressedLepton_HH
#define RIVET_DressedLepton_HH

#include "Rivet/Particle.hh"
#include "Rivet/Config/RivetCommon.hh"

namespace Rivet {


  /// Accepted classes of lepton origin
  enum class LeptonOrigin { PROMPT=1, NODECAY=1, ALL };

  /// Reconstruction/dressing mode for leptons
  enum class LeptonReco { ALL=0, ALL_BARE=0, ALL_DRESSED, PROMPT_BARE, PROMPT_DRESSED };

  /// The approach taken to photon dressing of leptons
  enum class DressingType { DR=0, CONE=0, CLUSTER=1, AKT=1 };

  /// Accepted classes of lepton origin  
  enum class PhotonOrigin { NONE=0, PROMPT=1, NODECAY=1, ALL };

  /// @brief Whether to add photons as resonance constituents
  ///
  /// @deprecated Immediately to be removed, as it doesn't do what it
  /// says, and is a dreadful idea better achieved through modern
  /// Particle::constituents()
  enum class PhotonsAsConstituents { NO, YES };



  /// @brief A charged lepton meta-particle created by clustering photons close to the bare lepton
  ///
  /// @deprecated Prefer to use Particle.constituents()
  class DressedLepton : public Particle {
  public:

    /// Copy constructor (from Particle)
    DressedLepton(const Particle& dlepton);

    /// @brief Components constructor
    ///
    /// @note This is not a copy constructor, hence the explicit second argument even if empty
    DressedLepton(const Particle& lepton, const Particles& photons, bool momsum=true);

    /// @brief Add a photon to the dressed lepton
    ///
    /// @todo Deprecate and override add/setConstituents instead?
    void addPhoton(const Particle& p, bool momsum=true);

    /// Retrieve the bare lepton
    const Particle& bareLepton() const;

    /// Retrieve the clustered photons
    const Particles photons() const { return slice(constituents(), 1); }

  };


  /// Alias for a list of dressed leptons, cf. Particles and Jets
  using DressedLeptons = vector<DressedLepton>;
  
  
}

#endif
