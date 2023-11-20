// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {
 

  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// A convenience method for finding l+l- pair resonances, including photon clustering.
  class ZFinder : public ParticleFinder {
  public:

    /// @name Constructors
    /// @{

    /// @brief Constructor taking cuts object
    ///
    /// @param inputfs Input final state
    /// @param cuts  Lepton cuts
    /// @param pid  Type of the leptons
    /// @param minmass,maxmass  Dilepton mass window
    /// @param chLeptons  The type of charged leptons considered
    /// @param dRmax  Maximum dR of photons around leptons to take into account
    ///  for Z reconstruction.
    /// @param clusterPhotons  Whether such photons are supposed to be
    ///  clustered to the lepton objects and thus Z mom
    /// @param masstarget  The expected (transverse) mass value, if resolving ambiguities
    ZFinder(const FinalState& inputfs,
		const Cut& cuts,
		PdgId pid,
		double minmass, double maxmass,
		double dRmax=0.1,
		LeptonOrigin chLeptons=LeptonOrigin::PROMPT,
	        PhotonOrigin clusterPhotons=PhotonOrigin::NODECAY,
	        double masstarget=91.2*GeV);


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(ZFinder);

    /// @}


    /// @brief Access to the found bosons
    ///
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& bosons() const { return particles(); }
    /// Access to the found boson (assuming it exists).
    const Particle& boson() const { return bosons().front(); }


    /// @brief Access to the constituent clustered leptons (and photons)
    ///
    /// The positive charge constituent is first in the list (if not empty), and
    /// the negative one second. Any included photons follow.
    const Particles& constituents() const;

    /// @brief Access to the constituent clustered leptons (and photons)
    ///
    /// As for constituents(), but returning 0 or 2 leptons only.
    Particles leptons() const {
      return head(constituents(), 2);
    }

    
    /// Access to the event-particles other than the Z leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() { _theParticles.clear(); }

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


  protected:

    /// Mass cuts to apply to clustered leptons (cf. InvMassFinalState)
    double _minmass, _maxmass, _masstarget;

    /// Lepton flavour
    PdgId _pid;

  };


  // /// Legacy name
  // using ZFinder = DilepFinder;


}

#endif
