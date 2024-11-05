// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Jet.fhh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetFastJet.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include <numeric>

namespace Rivet {


  /// @brief Specialised vector of Jet objects.
  ///
  /// A specialised version of vector<Jet> which is able to implicitly and
  /// explicitly convert to a vector of FourMomentum or PseudoJet.
  class Jets : public std::vector<Jet> {
  public:
    using base = std::vector<Jet>; //< using-declarations don't like template syntax
    using base::base; //< import base-class constructors
    Jets();
    Jets(const std::vector<Jet>& vjs);
    FourMomenta moms() const;
    PseudoJets pseudojets() const;
    operator FourMomenta () const { return moms(); }
    operator PseudoJets () const { return pseudojets(); }
    Jets& operator += (const Jet& j);
    Jets& operator += (const Jets& js);
  };

  Jets operator + (const Jets& a, const Jets& b);


  //////////////////////


  /// @brief Representation of a clustered jet of particles.
  class Jet : public ParticleBase {
  public:

    /// @name Constructors
    /// @{

    /// Constructor from a FastJet PseudoJet, with optional full particle constituents information.
    Jet(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pj, particles, tags);
    }

    /// Set the jet data, with optional full particle information.
    Jet(const FourMomentum& pjet, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pjet, particles, tags);
    }

    /// Default constructor -- only for STL storability
    Jet() { clear(); }

    /// @}


    /// @name Access jet constituents
    /// @{

    /// Number of particles in this jet.
    size_t size() const { return _particles.size(); }

    /// Get the particles in this jet.
    Particles& particles() { return _particles; }
    /// Get the particles in this jet (const version)
    const Particles& particles() const { return _particles; }
    /// Get the particles in this jet which pass a cut (const)
    const Particles particles(const Cut& c) const { return select(_particles, c); }
    /// Get the particles in this jet which pass a filtering functor (const)
    const Particles particles(const ParticleSelector& s) const { return select(_particles, s); }

    /// Get the particles in this jet (FastJet-like alias)
    Particles& constituents() { return particles(); }
    /// Get the particles in this jet (FastJet-like alias, const version)
    const Particles& constituents() const { return particles(); }
    /// Get the particles in this jet which pass a cut (FastJet-like alias, const)
    const Particles constituents(const Cut& c) const { return particles(c); }
    /// Get the particles in this jet which pass a filtering functor (FastJet-like alias, const)
    const Particles constituents(const ParticleSelector& s) const { return particles(s); }

    /// Check whether this jet contains a particular particle.
    bool containsParticle(const Particle& particle) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const Particle& particle) const { return containsParticle(particle); }

    /// Check whether this jet contains a certain particle type.
    bool containsParticleId(PdgId pid) const;
    /// Nicer alias for containsParticleId
    bool containsPID(PdgId pid) const { return containsParticleId(pid); }

    /// Check whether this jet contains at least one of certain particle types.
    bool containsParticleId(const vector<PdgId>& pids) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const vector<PdgId>& pids) const { return containsParticleId(pids); }

    /// @}


    /// @name Tagging
    ///
    /// @note General sources of tag particles are planned. The default jet finding
    /// adds b-hadron, c-hadron, and tau tags by ghost association.
    /// @{

    /// @brief Particles which have been tag-matched to this jet
    Particles& tags() { return _tags; }
    /// @brief Particles which have been tag-matched to this jet (const version)
    const Particles& tags() const { return _tags; }
    /// @brief Particles which have been tag-matched to this jet _and_ pass a selector function or dR requirement
    ///
    /// @note Note the less efficient return by value, due to the filtering.
    Particles tags(const ParticleSelector& f, double dRmax=-1) const {
      const Particles filttags = select(tags(), f);
      // Return, via dR filter if requested
      return dRmax < 0 ? filttags : select(filttags, deltaRLess(this->mom(), dRmax));
    }
    /// @brief Particles which have been tag-matched to this jet _and_ pass a Cut or dR requirement
    ///
    /// @note Note the less efficient return by value, due to the cut-pass filtering.
    Particles tags(const Cut& c, double dRmax=-1) const;


    /// @brief Get the b particles tag-matched to this jet
    ///
    /// The default jet finding adds b-hadron tags by ghost association.
    ///
    /// This version has an optional tag-particle Cut requirement
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles bTags(const Cut& c=Cuts::open(), double dRmax=-1) const;

    /// @brief Get the b particles tag-matched to this jet (with optional selector function and dR restriction)
    ///
    /// The default jet finding adds b-hadron tags by ghost association.
    ///
    /// This version has an optional tag-particle selector function
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles bTags(const ParticleSelector& f, double dRmax=-1) const { return select(bTags(), f); }

    /// Does this jet have at least one b-tag? (with optional Cut and dR restriction)
    bool bTagged(const Cut& c=Cuts::open(), double dRmax=-1) const { return !bTags(c).empty(); }

    /// Does this jet have at least one b-tag? (with optional selector function and dR restriction)
    bool bTagged(const ParticleSelector& f, double dRmax=-1) const { return !bTags(f).empty(); }


    /// @brief Get the c (and not b) particles tag-matched to this jet
    ///
    /// The default jet finding adds c-hadron tags by ghost association.
    ///
    /// This version has an optional tag-particle Cut requirement
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles cTags(const Cut& c=Cuts::open(), double dRmax=-1) const;

    /// @brief Get the c (and not b) particles which have been tag-matched to this jet
    ///
    /// The default jet finding adds c-hadron tags by ghost association.
    ///
    /// This version has an optional tag-particle selector function
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles cTags(const ParticleSelector& f, double dRmax=-1) const { return select(cTags(), f); }

    /// Does this jet have at least one c-tag? (with optional Cut and dR restriction)
    bool cTagged(const Cut& c=Cuts::open(), double dRmax=-1) const { return !cTags(c).empty(); }

    /// Does this jet have at least one c-tag? (with optional selector function and dR restriction)
    bool cTagged(const ParticleSelector& f, double dRmax=-1) const { return !cTags(f).empty(); }


    /// @brief Get the tau particles tag-matched to this jet
    ///
    /// The default jet finding adds tau tags by ghost association.
    ///
    /// This version has an optional tag-particle Cut requirement
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles tauTags(const Cut& c=Cuts::open(), double dRmax=-1) const;

    /// @brief Get the tau particles tag-matched to this jet
    ///
    /// The default jet finding adds tau tags by ghost association.
    ///
    /// This version has an optional tag-particle selector function
    /// and restriction to a tighter dR cone around the jet centroid.
    Particles tauTags(const ParticleSelector& f, double dRmax=-1) const { return select(tauTags(), f); }

    /// Does this jet have at least one tau-tag (with optional Cut and dR restriction)
    bool tauTagged(const Cut& c=Cuts::open(), double dRmax=-1) const { return !tauTags(c).empty(); }

    /// Does this jet have at least one tau-tag (with optional selector function and dR restriction)
    bool tauTagged(const ParticleSelector& f, double dRmax=-1) const { return !tauTags(f).empty(); }


    /// @todo Extend to arbitrary tagging types

    /// @}


    /// @name Effective jet 4-vector properties
    /// @{

    /// Get equivalent single momentum four-vector.
    const FourMomentum& momentum() const { return _momentum; }

    /// Apply an active Lorentz transform to this jet
    /// @note The Rivet jet momentum, constituent particles, and tag particles will be modified.
    /// @warning The FastJet cluster sequence and pseudojets will not be modified: don't use them after transformation!
    Jet& transformBy(const LorentzTransform& lt);

    /// Get the total energy of this jet.
    double totalEnergy() const { return momentum().E(); }

    /// Get the energy carried in this jet by neutral particles.
    double neutralEnergy() const;

    /// Get the energy carried in this jet by hadrons.
    double hadronicEnergy() const;

    /// @}


    /// @name Interaction with FastJet
    /// @{

    /// Access the internal FastJet3 PseudoJet (as a const reference)
    const fastjet::PseudoJet& pseudojet() const { return _pseudojet; }

    /// Cast operator to FastJet3 PseudoJet (as a const reference)
    operator const fastjet::PseudoJet& () const { return pseudojet(); }

    /// @}


    /// @name Set the jet constituents and properties
    /// @{

    /// @brief Set the jet data from a FastJet PseudoJet, with optional particle constituents and tags lists.
    ///
    /// @note The particles() list will be extracted from PseudoJet constituents
    /// by default, making use of an attached user info if one is found.
    Jet& setState(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles());

    /// Set all the jet data, with optional full particle constituent and tag information.
    Jet& setState(const FourMomentum& mom, const Particles& particles, const Particles& tags=Particles());

    /// @brief Set the particles collection with full particle information.
    ///
    /// If set, this overrides particle info extracted from the PseudoJet
    Jet& setParticles(const Particles& particles);
    Jet& setConstituents(const Particles& particles) { return setParticles(particles); }

    /// Reset this jet as empty.
    Jet& clear();

    /// @}


  private:

    /// FJ3 PseudoJet member to unify PseudoJet and Jet
    fastjet::PseudoJet _pseudojet;

    /// Full constituent particle information. (Filled from PseudoJet if possible.)
    /// @todo Make these mutable or similar? Add a flag to force a cache rebuild?
    Particles _particles;

    /// Particles used to tag this jet (can be anything, but c and b hadrons are the most common)
    Particles _tags;

    /// Effective jet 4-vector (just for caching)
    mutable FourMomentum _momentum;

    // /// Provide but hide the equality operators, to avoid implicit comparison via fastjet::PseudoJet
    // bool operator == (const Jet&) const;
    // bool operator != (const Jet&) const;

  };


  /// @name String representation and streaming support
  /// @{

  /// Allow a Jet to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const Jet& j);

  /// @}


}


#include "Rivet/Tools/JetUtils.hh"

#endif
