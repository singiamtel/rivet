// -*- C++ -*-
#ifndef RIVET_SmearedMET_HH
#define RIVET_SmearedMET_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/METFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing missing (transverse) energy/momentum with detector resolutions
  class SmearedMET : public METFinder {
  public:

    /// @name Constructors etc.
    /// @{

    /// @brief Constructor from a MissingMomentum projection and a smearing-params function
    ///
    /// Smearing-params function maps a 3-vector MET and scalar SET to a new METSmearParams struct for sampling
    template<typename SMEARPARAMSFN, typename std::enable_if_t<is_same_v<invoke_result_t<SMEARPARAMSFN, Vector3, double>, METSmearParams>, int> = 0>
    SmearedMET(const MissingMomentum& mm, const SMEARPARAMSFN& metSmearParamsFn)
      : _metSmearParamsFn(metSmearParamsFn), _metSmearFn(nullptr)
        // _metSmearFn([&](const Vector3& met, double set) -> Vector3 {
        //   const METSmearParams msps = metSmearParamsFn(met, set);
        //   return MET_SMEAR_NORM(msps);
        // })
    {
      setName("SmearedMET");
      declare(mm, "TruthMET");
    }

    /// @brief Constructor from a Cut (on the particles used to determine missing momentum) and a smearing-params function
    ///
    /// Smearing-params function maps a 3-vector MET and scalar SET to a new METSmearParams struct for sampling
    template<typename SMEARPARAMSFN, typename std::enable_if_t<is_same_v<invoke_result_t<SMEARPARAMSFN, Vector3, double>, METSmearParams>, int> = 0>
    SmearedMET(const SMEARPARAMSFN& metSmearParamsFn, const Cut& cut=Cuts::OPEN)
      : SmearedMET(MissingMomentum(cut), metSmearParamsFn)
    {  }

    /// @brief Constructor from a MissingMomentum projection and a smearing function
    ///
    /// Smearing function maps a 3-vector MET and scalar SET to a new MET 3-vector: f(V3, double) -> V3
    template<typename SMEARFN, typename std::enable_if_t<is_same_v<invoke_result_t<SMEARFN, Vector3, double>, Vector3>, int> = 0>
    SmearedMET(const MissingMomentum& mm, const SMEARFN& metSmearFn)
       : _metSmearParamsFn(nullptr), _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      declare(mm, "TruthMET");
    }

    /// @brief Constructor from a Cut (on the particles used to determine missing momentum) and a smearing function
    ///
    /// Smearing function maps a 3-vector MET and scalar SET to a new MET 3-vector: f(V3, double) -> V3
    template<typename SMEARFN, typename std::enable_if_t<is_same_v<invoke_result_t<SMEARFN, Vector3, double>, Vector3>, int> = 0>
    SmearedMET(const SMEARFN& metSmearFn, const Cut& cut=Cuts::OPEN)
     : SmearedMET(MissingMomentum(cut), metSmearFn)
    {  }

    /// @brief Constructor from a MissingMomentum projection and a pair of smearing-params and smearing functions
    ///
    /// Smearing-params function maps a 3-vector MET and scalar SET to a new METSmearParams struct for sampling
    ///
    /// Smearing function maps a 3-vector MET and scalar SET to a new MET 3-vector: f(V3, double) -> V3
    template <typename SMEARPARAMSFN, typename SMEARFN,
              typename std::enable_if_t<is_same_v<invoke_result_t<SMEARFN, Vector3, double>, Vector3> &&
                                        is_same_v<invoke_result_t<SMEARPARAMSFN, Vector3, double>, METSmearParams>, int> = 0>
    SmearedMET(const MissingMomentum& mm, const SMEARPARAMSFN& metSmearParamsFn, const SMEARFN& metSmearFn)
      : _metSmearParamsFn(metSmearParamsFn), _metSmearFn(metSmearFn)
    {
      setName("SmearedMET");
      declare(mm, "TruthMET");
    }

    /// @brief Constructor from a Cut (on the particles used to determine missing momentum) and a pair of smearing (params) functions
    ///
    /// Smearing-params function maps a 3-vector MET and scalar SET to a new METSmearParams struct for sampling
    ///
    /// Smearing function maps a 3-vector MET and scalar SET to a new MET 3-vector: f(V3, double) -> V3
    template <typename SMEARPARAMSFN, typename SMEARFN,
              typename std::enable_if_t<is_same_v<invoke_result_t<SMEARFN, Vector3, double>, Vector3> &&
                                        is_same_v<invoke_result_t<SMEARPARAMSFN, Vector3, double>, METSmearParams>, int> = 0>
    SmearedMET(const SMEARPARAMSFN& metSmearParamsFn, const SMEARFN& metSmearFn, const Cut& cut=Cuts::OPEN)
      : SmearedMET(MissingMomentum(cut), metSmearParamsFn, metSmearFn)
    {  }


    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(SmearedMET);

    /// @}

    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// Compare to another SmearedMET
    CmpState compare(const Projection& p) const {
      // const SmearedMET& other = dynamic_cast<const SmearedMET&>(p);
      // if (get_address(_metSmearParamsFn) == 0) return cmp((size_t)this, (size_t)&p);
      // if (get_address(_metSmearFn) == 0) return cmp((size_t)this, (size_t)&p);
      // MSG_TRACE("Smear hashes (params) = " << get_address(_metSmearParamsFn) << "," << get_address(other._metSmearParamsFn));
      // MSG_TRACE("Smear hashes (smear) = " << get_address(_metSmearFn) << "," << get_address(other._metSmearFn));
      // return mkPCmp(other, "TruthMET") ||
      //   cmp(get_address(_metSmearParamsFn), get_address(other._metSmearParamsFn));
      //   cmp(get_address(_metSmearFn), get_address(other._metSmearFn));
      return CmpState::UNDEF;
    }


    /// Perform the MET finding & smearing calculation
    void project(const Event& e) {
      const METFinder& mm = apply<MissingMomentum>(e, "TruthMET");
      _set = mm.scalarEt();
      _vet = mm.vectorEt();
      if (_metSmearFn) {
        _vet = _metSmearFn(_vet, mm.scalarEt()); //< custom smearing
      } else if (_metSmearParamsFn) {
        const METSmearParams msps = _metSmearParamsFn(_vet, mm.scalarEt()); //< custom smear params
        _vet = MET_SMEAR_NORM(msps); //< normal-distribution smearing from custom params
      } else {
        throw SmearError("Attempt to smear MET with neither smearing function nor smearing-params function set");
      }
    }


    /// @name Transverse-momentum functions
    ///
    /// @note This may be what you want, even if the paper calls it "missing Et"!
    /// @{

    /// The vector-summed visible transverse momentum in the event, as a 3-vector with z=0
    const Vector3& vectorPt() const { return vectorEt(); }
    /// The scalar-summed visible transverse momentum in the event
    double scalarPt() const { return scalarEt(); }

    /// @}


    /// @name Transverse-energy functions
    ///
    /// @warning Despite the common names "MET" and "SET", what's often meant is the pT functions above!
    /// @{

    /// The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing ET vector.
    const Vector3& vectorEt() const { return _vet; }

    /// The scalar-summed visible transverse energy in the event, as a 3-vector with z=0
    double scalarEt() const { return _set; }

    /// Obtain an approximation to the MET resolution for this event
    double missingEtResolution() const {
      if (!_metSmearParamsFn)
        throw UserError("Trying to compute MET significance without a registered significance function");
      METSmearParams msps = _metSmearParamsFn(vectorEt(), scalarEt());
      //return max(msps.eResolution, msps.pResolution);
      return msps.pResolution;
    }

    /// Obtain an approximation to the MET significant (value/resolution) for this event
    double missingEtSignf() const {
      return missingEt() / missingEtResolution();
    }

    /// @}


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() {  }


  protected:

    Vector3 _vet;
    double _set;

    /// Stored smearing-params function
    METSmearParamsFn _metSmearParamsFn;

    /// Stored smearing function
    METSmearFn _metSmearFn;

  };


}

#endif
