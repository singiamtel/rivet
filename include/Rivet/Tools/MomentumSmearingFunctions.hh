// -*- C++ -*-
#ifndef RIVET_MomentumSmearingFunctions_HH
#define RIVET_MomentumSmearingFunctions_HH

#include "Rivet/Math/Vector4.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @ingroup smearing
  /// @{

  /// @defgroup smearing_mom Generic 4-momentum filtering, efficiency and smearing utils
  /// @{

  /// @brief Struct for holding momentum-smearing parameters
  ///
  /// Contains a nominal 4-momentum (maybe shifted in energy or angle)
  /// and absolute resolutions in momentum, energy, phi, and eta, with
  /// which to smear the nominal.
  struct P4SmearParams {
    P4SmearParams(const FourMomentum& p4nom, double pres=0.0, double eres=0.0, double phires=0.0, double etares=0.0)
      : p4Nominal(p4nom), pResolution(pres), eResolution(pres),
        phiResolution(phires), etaResolution(etares) {}
    FourMomentum p4Nominal;
    double pResolution{0.0}, eResolution{0.0}, phiResolution{0.0}, etaResolution{0.0};
  };

  /// Typedef for FourMomentum smearing functions/functors
  typedef std::function<FourMomentum(const FourMomentum&)> P4SmearFn;

  /// Typedef for FourMomentum efficiency functions/functors
  typedef std::function<double(const FourMomentum&)> P4EffFn;


  /// Take a FourMomentum and return 0
  inline double P4_EFF_ZERO(const FourMomentum& ) { return 0; }

  /// Take a FourMomentum and return 1
  inline double P4_EFF_ONE(const FourMomentum& ) { return 1; }

  /// Take a FourMomentum and return a constant number
  struct P4_EFF_CONST {
    P4_EFF_CONST(double x) : _x(x) {}
    double operator () (const FourMomentum& )  const { return _x; }
    double _x;
  };


  /// Take a FourMomentum and return it unmodified
  inline FourMomentum P4_SMEAR_IDENTITY(const FourMomentum& p) { return p; }
  /// Alias for P4_SMEAR_IDENTITY
  inline FourMomentum P4_SMEAR_PERFECT(const FourMomentum& p) { return p; }

  /// Smear a FourMomentum's energy using a Gaussian of absolute width @a resolution
  /// @todo Also make jet versions that update/smear constituents?
  inline FourMomentum P4_SMEAR_E_GAUSS(const FourMomentum& p, double resolution) {
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    const double smeared_E = max(randnorm(p.E(), resolution), mass); //< can't let the energy go below the mass!
    return FourMomentum::mkEtaPhiME(p.eta(), p.phi(), mass, smeared_E);
  }

  /// Smear a FourMomentum's transverse momentum using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_PT_GAUSS(const FourMomentum& p, double resolution) {
    const double smeared_pt = max(randnorm(p.pT(), resolution), 0.);
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), mass, smeared_pt);
  }

  /// Smear a FourMomentum's mass using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_MASS_GAUSS(const FourMomentum& p, double resolution) {
    const double smeared_mass = max(randnorm(p.mass(), resolution), 0.);
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), smeared_mass, p.pT());
  }

  /// @}



  /// @defgroup smearing_mom Generic 3-momentum filtering, efficiency and smearing utils
  /// @{

  /// Take a Vector3 and return 0
  inline double P3_EFF_ZERO(const Vector3&) { return 0; }

  /// Take a Vector3 and return 1
  inline double P3_EFF_ONE(const Vector3&) { return 1; }

  /// Take a Vector3 and return a constant number
  struct P3_EFF_CONST {
    P3_EFF_CONST(double x) : _x(x) {}
    double operator () (const Vector3& )  const { return _x; }
    double _x;
  };


  /// Take a Vector3 and return it unmodified
  inline Vector3 P3_SMEAR_IDENTITY(const Vector3& p) { return p; }
  /// Alias for P3_SMEAR_IDENTITY
  inline Vector3 P3_SMEAR_PERFECT(const Vector3& p) { return p; }

  /// Smear a Vector3's length using a Gaussian of absolute width @a resolution
  inline Vector3 P3_SMEAR_LEN_GAUSS(const Vector3& p, double resolution) {
    const double smeared_mod = max(randnorm(p.mod(), resolution), 0.); //< can't let the energy go below the mass!
    return smeared_mod * p.unit();
  }

  /// @}



  /// @defgroup smearing_met Generic missing-ET smearing functions
  /// @{

  /// @brief Struct for holding MET-smearing parameters
  ///
  /// Contains a nominal 3-momentum (maybe shifted in energy or angle)
  /// and absolute resolutions in momentum, energy, phi, and eta, with
  /// which to Gaussian-smear the nominal.
  struct METSmearParams {
    METSmearParams(const Vector3& vmetnom, double pres=0.0, double phires=0.0)
      : vmetNominal(vmetnom), pResolution(pres), phiResolution(phires) {}
    Vector3 vmetNominal;
    double pResolution{0.0}, phiResolution{0.0};
  };

  /// Typedef for MET smearing-parameter functions/functors (given a MET vector and scalar sum(ET))
  ///
  /// Returns a nominal MET 3-vector (maybe shifted in energy or angle)
  /// and a resolution (in GeV) with which to Gaussian-smear its momentum.
  ///
  /// @todo Allow MET calculation to access the whole Event?
  typedef function<METSmearParams(const Vector3&, double)> METSmearParamsFn;

  /// Typedef for MET smearing functions/functors (given a MET vector and scalar sum(ET))
  ///
  /// @todo Allow MET calculation to access the whole Event?
  typedef function<Vector3(const Vector3&, double)> METSmearFn;


  /// @brief Smear a nominal vector magnitude by Gaussian with the given absolute resolutions
  ///
  /// @note For momentum resolutions comparable to or larger than the
  /// value they are smearing, a modulus of the resulting vaue is
  /// taken, i.e. "reflecting" negative tails back up to positive
  /// rather than truncating at zero (or -- horrors -- accidentally
  /// flipping the vector direction!) This is not obviously _better_
  /// than truncating at zero, and ideally we'd avoid this with some
  /// physical positive-definite sampling... but until then let's try
  /// to make sure the assuming-Gaussian resolution functions don't
  /// have to work in very asymmetric regimes where a single sigma
  /// isn't enough.
  ///
  /// @todo Improve on "reflect in zero" sampling behaviour where res >~ mu.
  inline Vector3 MET_SMEAR_NORM(const METSmearParams& msps) {
    const Vector3& vmet = msps.vmetNominal;
    // MET magnitude smear
    // const double metsmear = max(randnorm(vmet.mod(), resolution), 0.); //< make peak at 0
    const double metsmear = fabs(randnorm(vmet.mod(), msps.pResolution)); //< "reflect" at 0
    // MET phi-angle smearing matrix
    const Matrix3 phismear = Matrix3::mkZRotation(randnorm(0.0, msps.phiResolution));
    // Apply smearings
    return metsmear * phismear * vmet.unit();
  }


  /// @brief Identity resolution is 0 (perfect delta function, no offset)
  ///
  /// @note This simple factorization doesn't allow for expression of
  /// offsets or angular smearing... revisit if needed
  inline METSmearParams MET_SMEARPARAMS_IDENTITY(const Vector3& met, double) {
    return METSmearParams(met, 0.0, 0.0);
  }

  /// @brief Identity MET smearing just returns the input
  inline Vector3 MET_SMEAR_IDENTITY(const Vector3& met, double) {
    return met;
  }

  /// @}

  /// @}

}

#endif
