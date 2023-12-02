// -*- C++ -*-
#ifndef RIVET_DileptonFinder_HH
#define RIVET_DileptonFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// A convenience method for finding l+l- pair resonances, including photon clustering.
  class DileptonFinder : public FinalState {
  public:

    /// @name Constructors
    /// @{

    /// Modern constructor with a default FS, explicit mass target and dressing dR up-front (PID and LL via cuts)
    DileptonFinder(double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   const Cut& llcuts=Cuts::OPEN,
		   LeptonOrigin leptonOrigin=LeptonOrigin::PROMPT,
		   PhotonOrigin photonOrigin=PhotonOrigin::NODECAY,
                   TauDecaysAs tauDecays=TauDecaysAs::PROMPT,
                   MuDecaysAs muDecays=MuDecaysAs::PROMPT,
                   DressingType dressing=DressingType::CONE);


    /// Modern constructor with an explicit FS, explicit mass target and dressing dR up-front (PID and LL via cuts)
    DileptonFinder(const FinalState& allfs,
                   double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   const Cut& llcuts=Cuts::OPEN,
                   DressingType dressing=DressingType::CONE);


    /// Modern constructor with two explicit FSes, explicit mass target and dressing dR up-front (PID and LL via cuts)
    DileptonFinder(const FinalState& leptonfs,
                   const FinalState& photonfs,
                   double masstarget,
                   double dRmax,
                   const Cut& lcuts=Cuts::OPEN,
                   const Cut& llcuts=Cuts::OPEN,
                   DressingType dressing=DressingType::CONE);


    /// @brief Old constructor taking min and max ll masses
    ///
    /// @deprecated Avoid!
    DileptonFinder(const FinalState& allfs,
		   const Cut& lcuts,
		   double masstarget,
                   PdgId pid,
		   double minmass, double maxmass,
		   double dRdress)
      : DileptonFinder(allfs, masstarget, dRdress,
                       lcuts && Cuts::abspid == pid,
                       Cuts::mass >= minmass && Cuts::mass < maxmass)
    {   }


    /// Ok constructor with an explicit mass target and dressing dR up-front
    ///
    /// @deprecated Avoid!
    DileptonFinder(const FinalState& allfs,
		   PdgId pid,
		   double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   double minmass=-DBL_MAX, double maxmass=DBL_MAX)
      : DileptonFinder(allfs, masstarget, dRdress,
                       lcuts && Cuts::abspid == pid,
                       Cuts::mass >= minmass && Cuts::mass < maxmass)
    {   }


    /// Ok constructor with a default FS, explicit mass target and dressing dR up-front
    ///
    /// @deprecated Avoid!
    DileptonFinder(PdgId pid,
		   double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   double minmass=-DBL_MAX, double maxmass=DBL_MAX,
		   LeptonOrigin whichleptons=LeptonOrigin::PROMPT,
		   PhotonOrigin whichphotons=PhotonOrigin::NODECAY,
                   TauDecaysAs tauDecays=TauDecaysAs::PROMPT,
                   MuDecaysAs muDecays=MuDecaysAs::PROMPT,
                   DressingType dressing=DressingType::CONE)
      : DileptonFinder(masstarget, dRdress,
                       lcuts && Cuts::abspid == pid,
                       Cuts::mass >= minmass && Cuts::mass < maxmass,
                       whichleptons, whichphotons, tauDecays, muDecays, dressing)
    {   }


    /// Ok constructor with a default FS, explicit mass target and dressing dR up-front
    ///
    /// @deprecated Avoid!
    DileptonFinder(double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   double minmass=-DBL_MAX, double maxmass=DBL_MAX,
		   LeptonOrigin whichleptons=LeptonOrigin::PROMPT,
		   PhotonOrigin whichphotons=PhotonOrigin::NODECAY,
                   TauDecaysAs tauDecays=TauDecaysAs::PROMPT,
                   MuDecaysAs muDecays=MuDecaysAs::PROMPT,
                   DressingType dressing=DressingType::CONE)
      : DileptonFinder(masstarget, dRdress, lcuts,
                       Cuts::mass >= minmass && Cuts::mass < maxmass,
                       whichleptons, whichphotons, tauDecays, muDecays, dressing)
    {   }


    /// Ok constructor with an explicit FS, explicit mass target and dressing dR up-front (PID via cuts)
    ///
    /// @deprecated Avoid!
    DileptonFinder(const FinalState& allfs,
                   double masstarget,
		   double dRdress,
		   const Cut& lcuts=Cuts::OPEN,
		   double minmass=-DBL_MAX, double maxmass=DBL_MAX)
      : DileptonFinder(allfs, masstarget, dRdress, lcuts,
                       Cuts::mass >= minmass && Cuts::mass < maxmass)
    {   }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(DileptonFinder);

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
    ///
    /// The positive charge constituent is first in the list (if not empty), and
    /// the negative one second. Any included photons follow.
    Particles leptons() const {  return constituents(); }


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

    /// Mass target for identifying the best fit
    double _masstarget;

  };


  /// Aliases
  using ZFinder = DileptonFinder;
  using LLFinder = DileptonFinder;


}

#endif
