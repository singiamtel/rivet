// -*- C++ -*-
#ifndef RIVET_LeptonFinder_HH
#define RIVET_LeptonFinder_HH

#include "Rivet/DressedLepton.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Config/RivetCommon.hh"

namespace Rivet {
  

  /// @brief Reconstruct leptons, generally including "dressing" with clustered photons
  ///
  /// The clustering is done by a delta(R) cone around each bare lepton or by
  /// jet clustering. In both modes, double counting is avoided: for the dR
  /// clustering, a photon is only added to the _closest_ bare lepton if it
  /// happens to be within the capture radius of more than one; for the jet
  /// clustering, only the bare lepton with the highest pT is retained if more
  /// than one is clustered into a jet.
  ///
  /// @note The particles() and dressedLeptons() methods both return the
  /// composite clustered-lepton objects, just with a few extra helper methods
  /// on the special DressedLepton type returned by the latter. The constituent
  /// bare leptons and photons are returned by rawParticles() (inherited from
  /// ParticleFinder)
  ///
  class LeptonFinder : public FinalState {
  public:

    /// @brief Constructor with a single input FinalState (used for both photons and bare leptons)
    ///
    /// Provide a single final state projection used to select the photons and
    /// bare leptons, a photon-clustering delta(R) cone size around each bare
    /// lepton, and an optional cut on the _dressed_ leptons (i.e. the momenta
    /// and PID after clustering).  The final arguments control whether
    /// non-prompt photons are to be included, and whether the matching of
    /// photons to leptons is to be done via dR matching to the bare lepton or
    /// by a jet clustering algorithm.  Set the clustering radius to 0 or
    /// negative to disable clustering.
    LeptonFinder(const FinalState& allfs,
		 double dRmax, const Cut& cut=Cuts::OPEN,
		 PhotonOrigin whichphotons=PhotonOrigin::NODECAY,
		 DressingType dressing=DressingType::CONE);

    /// @brief Constructor with default input FinalState
    ///
    /// DressedLepton construction from a default-constructed FinalState.
    /// Provide a photon-clustering delta(R) cone size around each bare lepton,
    /// and an optional cut on the _dressed_ leptons (i.e. the momenta and PID
    /// after clustering).  The final arguments control whether non-prompt
    /// photons are to be included, and whether the matching of photons to
    /// leptons is to be done via dR matching to the bare lepton or by a jet
    /// clustering algorithm.  Set the clustering radius to 0 or negative to
    /// disable clustering.
    LeptonFinder(double dRmax, const Cut& cut=Cuts::OPEN,
		 PhotonOrigin whichphotons=PhotonOrigin::NODECAY,
		 DressingType dressing=DressingType::CONE)
      : LeptonFinder(FinalState(), dRmax, cut, whichphotons, dressing)
    {   }

    /// @brief Constructor with distinct photon and lepton finders
    ///
    /// Provide final state projections used to select the photons and bare
    /// leptons, a clustering delta(R) cone size around each bare lepton, and an
    /// optional cut on the _dressed_ leptons (i.e. the momenta and PID after
    /// clustering.)  The final arguments control whether non-prompt photons are
    /// to be included, and whether the matching of photons to leptons is to be
    /// done via dR matching to the bare lepton or by a jet clustering
    /// algorithm.  Set the clustering radius to 0 or negative to disable
    /// clustering.
    ///
    /// @note Wish we had put the first two args the other way around... not reversible!
    ///
    /// @todo Convert the second arg to a general ParticleFinder rather than an FS, to
    /// allow clustering of unstables, e.g. taus via TauFinder when that becomes a PF.
    /// Complicated by the clustering version relying on MergedFinalState and FastJets'
    /// current restriction to FinalState inputs. Requires widespread redesign.
    LeptonFinder(const FinalState& photons, const FinalState& bareleptons,
		 double dRmax, const Cut& cut=Cuts::OPEN,
		 PhotonOrigin whichphotons=PhotonOrigin::NODECAY,
		 DressingType dressing=DressingType::CONE);


    /// Clone this projection
    RIVET_DEFAULT_PROJ_CLONE(LeptonFinder);

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// @brief Retrieve the dressed leptons
    ///
    /// @note Like particles() but with helper functions
    DressedLeptons dressedLeptons() const {
      DressedLeptons rtn;
      for (const Particle& p : particles(cmpMomByPt))
        rtn += DressedLepton(p);  //static_cast<const DressedLepton>(p);
      return rtn;
    }

    /// @brief Retrieve the dressed leptons ordered by supplied sorting functor
    ///
    /// @note Like particles() but with helper functions
    DressedLeptons dressedLeptons(const ParticleSorter& sorter) const {
      DressedLeptons rtn;
      for (const Particle& p : particles(sorter))
        rtn += DressedLepton(p);  //static_cast<const DressedLepton>(p);
      return rtn;
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  protected:

    /// Maximum cone radius to find photons in
    double _dRmax;

    /// Whether to include photons from hadron (particularly pi0) and hadronic tau decays
    bool _fromDecay;

    /// Whether to use a jet clustering algorithm rather than nearest-lepton association
    bool _useJetClustering;


  };



}

#endif
