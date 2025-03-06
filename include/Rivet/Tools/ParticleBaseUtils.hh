#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/Tools/Utils.hh"
#include "Rivet/ParticleBase.hh"

namespace Rivet {


  /// @defgroup particlebaseutils Functions for Particles and Jets
  /// @{

  /// @defgroup particlebasetutils_pb2bool ParticleBase classifier -> bool functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to any() or all() e.g. any(jets, DeltaRLess(electron, 0.4))
  /// @{

  /// std::function instantiation for functors taking a ParticleBase and returning a bool
  using ParticleBaseSelector = function<bool(const ParticleBase&)>;
  /// std::function instantiation for functors taking two ParticleBase and returning a bool
  using ParticleBaseSorter = function<bool(const ParticleBase&, const ParticleBase&)>;


  /// Base type for Particle -> bool functors
  struct BoolParticleBaseFunctor {
    virtual bool operator()(const ParticleBase& p) const = 0;
    virtual ~BoolParticleBaseFunctor() {}
  };


  /// Transverse momentum greater-than functor
  struct PtGtr : public BoolParticleBaseFunctor {
    PtGtr(double pt) : ptcut(pt) { }
    PtGtr(const FourMomentum& p) : ptcut(p.pT()) { }
    bool operator()(const ParticleBase& p) const { return p.pT() > ptcut; }
    double ptcut;
  };
  using pTGtr = PtGtr;
  using ptGtr = PtGtr;

  /// Transverse momentum less-than functor
  struct PtLess : public BoolParticleBaseFunctor {
    PtLess(const FourMomentum& p) : ptcut(p.pT()) { }
    PtLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };
  using pTLess = PtLess;
  using ptLess = PtLess;

  /// Transverse momentum in-range functor
  struct PtInRange : public BoolParticleBaseFunctor {
    PtInRange(pair<double, double> ptcuts) : ptcut(ptcuts) { }
    PtInRange(double ptlow, double pthigh) : PtInRange(make_pair(ptlow, pthigh)) { }
    PtInRange(const FourMomentum& p1, const FourMomentum& p2) : PtInRange(p1.pT(), p2.pT()) { }
    bool operator()(const ParticleBase& p) const { return p.pT() >= ptcut.first && p.pT() < ptcut.second; }
    pair<double,double> ptcut;
  };
  using pTInRange = PtInRange;
  using ptInRange = PtInRange;


  /// Pseudorapidity greater-than functor
  struct EtaGtr : public BoolParticleBaseFunctor {
    EtaGtr(double eta) : etacut(eta) { }
    EtaGtr(const FourMomentum& p) : etacut(p.eta()) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };
  using etaGtr = EtaGtr;

  /// Pseudorapidity less-than functor
  struct EtaLess : public BoolParticleBaseFunctor {
    EtaLess(double eta) : etacut(eta) { }
    EtaLess(const FourMomentum& p) : etacut(p.eta()) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };
  using etaLess = EtaLess;

  /// Pseudorapidity in-range functor
  struct EtaInRange : public BoolParticleBaseFunctor {
    EtaInRange(pair<double, double> etacuts) : etacut(etacuts) { }
    EtaInRange(double etalow, double etahigh) : EtaInRange(make_pair(etalow, etahigh)) { }
    EtaInRange(const FourMomentum& p1, const FourMomentum& p2) : EtaInRange(p1.eta(), p2.eta()) { }
    bool operator()(const ParticleBase& p) const { return p.eta() >= etacut.first && p.eta() < etacut.second; }
    pair<double,double> etacut;
  };
  using etaInRange = EtaInRange;


  /// Abs pseudorapidity greater-than functor
  struct AbsEtaGtr : public BoolParticleBaseFunctor {
    AbsEtaGtr(double abseta) : absetacut(abseta) { }
    AbsEtaGtr(const FourMomentum& p) : absetacut(p.abseta()) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };
  using absEtaGtr = AbsEtaGtr;
  using absetaGtr = AbsEtaGtr;

  /// Abs pseudorapidity momentum less-than functor
  struct AbsEtaLess : public BoolParticleBaseFunctor {
    AbsEtaLess(double abseta) : absetacut(abseta) { }
    AbsEtaLess(const FourMomentum& p) : absetacut(p.abseta()) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };
  using absEtaLess = AbsEtaLess;
  using absetaLess = AbsEtaLess;

  /// Abs pseudorapidity in-range functor
  struct AbsEtaInRange : public BoolParticleBaseFunctor {
    AbsEtaInRange(const pair<double, double>& absetacuts) : absetacut(absetacuts) { }
    AbsEtaInRange(double absetalow, double absetahigh) : AbsEtaInRange(make_pair(absetalow, absetahigh)) { }
    AbsEtaInRange(const FourMomentum& p1, const FourMomentum& p2) : AbsEtaInRange(p1.abseta(), p2.abseta()) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() >= absetacut.first && p.abseta() < absetacut.second; }
    pair<double,double> absetacut;
  };
  using absEtaInRange = AbsEtaInRange;
  using absetaInRange = AbsEtaInRange;


  /// Rapidity greater-than functor
  struct RapGtr : public BoolParticleBaseFunctor {
    RapGtr(double rap) : rapcut(rap) { }
    RapGtr(const FourMomentum& p) : rapcut(p.rap()) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };
  using rapGtr = RapGtr;

  /// Rapidity momentum less-than functor
  struct RapLess : public BoolParticleBaseFunctor {
    RapLess(double rap) : rapcut(rap) { }
    RapLess(const FourMomentum& p) : rapcut(p.rap()) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };
  using rapLess = RapLess;

  /// Rapidity in-range functor
  struct RapInRange : public BoolParticleBaseFunctor {
    RapInRange(const pair<double, double>& rapcuts) : rapcut(rapcuts) { }
    RapInRange(double raplow, double raphigh) : RapInRange(make_pair(raplow, raphigh)) { }
    RapInRange(const FourMomentum& p1, const FourMomentum& p2) : RapInRange(p1.rap(), p2.rap()) { }
    bool operator()(const ParticleBase& p) const { return p.rap() >= rapcut.first && p.rap() < rapcut.second; }
    pair<double,double> rapcut;
  };
  using rapInRange = RapInRange;


  /// Abs rapidity greater-than functor
  struct AbsRapGtr : public BoolParticleBaseFunctor {
    AbsRapGtr(double absrap) : absrapcut(absrap) { }
    AbsRapGtr(const FourMomentum& p) : absrapcut(p.absrap()) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };
  using absRapGtr = AbsRapGtr;
  using absrapGtr = AbsRapGtr;

  /// Abs rapidity momentum less-than functor
  struct AbsRapLess : public BoolParticleBaseFunctor {
    AbsRapLess(double absrap) : absrapcut(absrap) { }
    AbsRapLess(const FourMomentum& p) : absrapcut(p.absrap()) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };
  using absRapLess = AbsRapLess;
  using absrapLess = AbsRapLess;

  /// Abs rapidity in-range functor
  struct AbsRapInRange : public BoolParticleBaseFunctor {
    AbsRapInRange(const pair<double, double>& absrapcuts) : absrapcut(absrapcuts) { }
    AbsRapInRange(double absraplow, double absraphigh) : AbsRapInRange(make_pair(absraplow, absraphigh)) { }
    AbsRapInRange(const FourMomentum& p1, const FourMomentum& p2) : AbsRapInRange(p1.absrap(), p2.absrap()) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() >= absrapcut.first && p.absrap() < absrapcut.second; }
    pair<double,double> absrapcut;
  };
  using absRapInRange = AbsRapInRange;
  using absrapInRange = AbsRapInRange;


  /// @todo Define dR and dphi functors w.r.t. *multiple* ref vectors, with "all" or "any" semantics


  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) greater-than functor
  struct DeltaRGtr : public BoolParticleBaseFunctor {
    DeltaRGtr(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) > drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRGtr = DeltaRGtr;

  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) less-than functor
  struct DeltaRLess : public BoolParticleBaseFunctor {
    DeltaRLess(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) < drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRLess = DeltaRLess;

  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaRInRange : public BoolParticleBaseFunctor {
    DeltaRInRange(const ParticleBase& vec, const pair<double,double>& dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRInRange(const ParticleBase& vec, double drmin, double drmax, RapScheme scheme=PSEUDORAPIDITY)
      : DeltaRInRange(vec, make_pair(drmin, drmax), scheme) { }
    DeltaRInRange(const FourMomentum& vec, const pair<double,double>& dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRInRange(const FourMomentum& vec, double drmin, double drmax, RapScheme scheme=PSEUDORAPIDITY)
      : DeltaRInRange(vec, make_pair(drmin, drmax), scheme) { }
    DeltaRInRange(const Vector3& vec, const pair<double,double>& dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    DeltaRInRange(const Vector3& vec, double drmin, double drmax)
      : DeltaRInRange(vec, make_pair(drmin, drmax)) { }
    bool operator()(const ParticleBase& p) const {
      const double dR = deltaR(p, refvec, rapscheme);
      return dR >= drcut.first && dR < drcut.second;
    }
    FourMomentum refvec;
    pair<double,double> drcut;
    RapScheme rapscheme;
  };
  using deltaRInRange = DeltaRInRange;


  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaPhiGtr : public BoolParticleBaseFunctor {
    DeltaPhiGtr(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) > dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiGtr = DeltaPhiGtr;

  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaPhiLess : public BoolParticleBaseFunctor {
    DeltaPhiLess(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) < dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiLess = DeltaPhiLess;

  /// @f$ \Delta \phi @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaPhiInRange : public BoolParticleBaseFunctor {
    DeltaPhiInRange(const ParticleBase& vec, const pair<double,double>& dphi)
      : refvec(vec.mom()), dphicut(dphi) { }
    DeltaPhiInRange(const ParticleBase& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    DeltaPhiInRange(const FourMomentum& vec, const pair<double,double>& dphi)
      : refvec(vec), dphicut(dphi) { }
    DeltaPhiInRange(const FourMomentum& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    DeltaPhiInRange(const Vector3& vec, const pair<double,double>& dphi)
      : refvec(vec), dphicut(dphi) { }
    DeltaPhiInRange(const Vector3& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    bool operator()(const ParticleBase& p) const {
      const double dphi = deltaPhi(p, refvec);
      return dphi >= dphicut.first && dphi < dphicut.second;
    }
    Vector3 refvec;
    pair<double,double> dphicut;
  };
  using deltaPhiInRange = DeltaPhiInRange;


  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaEtaGtr : public BoolParticleBaseFunctor {
    DeltaEtaGtr(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) > detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaGtr = DeltaEtaGtr;

  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaEtaLess : public BoolParticleBaseFunctor {
    DeltaEtaLess(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) < detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaLess = DeltaEtaLess;

  /// @f$ \Delta \eta @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaEtaInRange : public BoolParticleBaseFunctor {
    DeltaEtaInRange(const ParticleBase& vec, const pair<double,double>& deta)
      : refvec(vec.mom()), detacut(deta) { }
    DeltaEtaInRange(const ParticleBase& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    DeltaEtaInRange(const FourMomentum& vec, const pair<double,double>& deta)
      : refvec(vec), detacut(deta) { }
    DeltaEtaInRange(const FourMomentum& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    DeltaEtaInRange(const Vector3& vec, const pair<double,double>& deta)
      : refvec(vec), detacut(deta) { }
    DeltaEtaInRange(const Vector3& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    bool operator()(const ParticleBase& p) const {
      const double deta = deltaEta(p, refvec);
      return deta >= detacut.first && deta < detacut.second;
    }
    Vector3 refvec;
    pair<double,double> detacut;
  };
  using deltaEtaInRange = DeltaEtaInRange;


  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaRapGtr : public BoolParticleBaseFunctor {
    DeltaRapGtr(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapGtr(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) > drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapGtr = DeltaRapGtr;

  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaRapLess : public BoolParticleBaseFunctor {
    DeltaRapLess(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapLess(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) < drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapLess = DeltaRapLess;

  /// @f$ \Delta y @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaRapInRange : public BoolParticleBaseFunctor {
    DeltaRapInRange(const ParticleBase& vec, const pair<double,double>& drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapInRange(const ParticleBase& vec, double drapmin, double drapmax)
      : DeltaRapInRange(vec, make_pair(drapmin, drapmax)) { }
    DeltaRapInRange(const FourMomentum& vec, const pair<double,double>& drap)
      : refvec(vec), drapcut(drap) { }
    DeltaRapInRange(const FourMomentum& vec, double drapmin, double drapmax)
      : DeltaRapInRange(vec, make_pair(drapmin, drapmax)) { }
    bool operator()(const ParticleBase& p) const {
      const double drap = deltaRap(p, refvec);
      return drap >= drapcut.first && drap < drapcut.second;
    }
    FourMomentum refvec;
    pair<double,double> drapcut;
  };
  using deltaRapInRange = DeltaRapInRange;

  /// @}


  /// @defgroup particlebaseutils_pb2dbl ParticleBase comparison -> double functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to transform()any(jets, DeltaRLess(electron, 0.4))
  /// @{

  /// Base type for Particle -> double functors
  struct DoubleParticleBaseFunctor {
    virtual double operator()(const ParticleBase& p) const = 0;
    virtual ~DoubleParticleBaseFunctor() {}
  };

  /// Calculator of @f$ \Delta R @f$ with respect to a given momentum
  struct DeltaRWRT : public DoubleParticleBaseFunctor {
    DeltaRWRT(const ParticleBase& pb, RapScheme scheme=PSEUDORAPIDITY) : p(pb.mom()), rapscheme(scheme) {}
    DeltaRWRT(const FourMomentum& p4, RapScheme scheme=PSEUDORAPIDITY) : p(p4), rapscheme(scheme) {}
    DeltaRWRT(const Vector3& p3) : p(p3.mod(), p3.x(), p3.y(), p3.z()), rapscheme(PSEUDORAPIDITY) {}
    double operator()(const ParticleBase& pb) const { return deltaR(p, pb, rapscheme); }
    double operator()(const FourMomentum& p4) const { return deltaR(p, p4, rapscheme); }
    double operator()(const Vector3& p3) const { return deltaR(p, p3); }
    const FourMomentum p;
    RapScheme rapscheme;
  };
  using deltaRWRT = DeltaRWRT;

  /// Calculator of @f$ \Delta \phi @f$ with respect to a given momentum
  struct DeltaPhiWRT : public DoubleParticleBaseFunctor {
    DeltaPhiWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaPhiWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaPhiWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaPhi(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaPhi(p, p4); }
    double operator()(const Vector3& p3) const { return deltaPhi(p, p3); }
    const Vector3 p;
  };
  using deltaPhiWRT = DeltaPhiWRT;

  /// Calculator of @f$ \Delta \eta @f$ with respect to a given momentum
  struct DeltaEtaWRT : public DoubleParticleBaseFunctor {
    DeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaEta(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaEta(p, p4); }
    double operator()(const Vector3& p3) const { return deltaEta(p, p3); }
    const Vector3 p;
  };
  using deltaEtaWRT = DeltaEtaWRT;

  /// Calculator of @f$ |\Delta \eta| @f$ with respect to a given momentum
  struct AbsDeltaEtaWRT : public DoubleParticleBaseFunctor {
    AbsDeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    AbsDeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    AbsDeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaEta(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaEta(p, p4)); }
    double operator()(const Vector3& p3) const { return fabs(deltaEta(p, p3)); }
    const Vector3 p;
  };
  using absDeltaEtaWRT = AbsDeltaEtaWRT;

  /// Calculator of @f$ \Delta y @f$ with respect to a given momentum
  struct DeltaRapWRT : public DoubleParticleBaseFunctor {
    DeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    DeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return deltaRap(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaRap(p, p4); }
    const FourMomentum p;
  };
  using deltaRapWRT = DeltaRapWRT;

  /// Calculator of @f$ |\Delta y| @f$ with respect to a given momentum
  struct AbsDeltaRapWRT : public DoubleParticleBaseFunctor {
    AbsDeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    AbsDeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaRap(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaRap(p, p4)); }
    const FourMomentum p;
  };
  using absDeltaRapWRT = AbsDeltaRapWRT;

  /// @}


  /// @defgroup particlebaseutils_uberfilt Next-level filtering
  /// @{

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void idiscardIfAny(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                            typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                        const typename PBCONTAINER2::value_type&)> fn) {
    for (const auto& pbcmp : tocompare) {
      idiscard(tofilter, [&](const typename PBCONTAINER1::value_type& pbfilt){ return fn(pbfilt, pbcmp); });
    }
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 discardIfAny(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                                   typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                               const typename PBCONTAINER2::value_type&)> fn) {
    PBCONTAINER1 tmp{tofilter};
    idiscardIfAny(tmp, tocompare, fn);
    return tmp;
  }


  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 selectIfAny(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                                  typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                              const typename PBCONTAINER2::value_type&)> fn) {
    PBCONTAINER1 selected;
    for (const auto& pbfilt : tofilter) {
      if (any(tocompare, [&](const typename PBCONTAINER2::value_type& pbcmp){ return fn(pbfilt, pbcmp); })) {
        selected += pbfilt;
      }
    }
    return selected;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void iselectIfAny(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                           typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                       const typename PBCONTAINER2::value_type&)> fn) {
    tofilter = selectIfAny(tofilter, tocompare, fn);
  }



  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 discardIfAll(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                                   typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                               const typename PBCONTAINER2::value_type&)> fn) {
    PBCONTAINER1 selected;
    for (const auto& pbfilt : tofilter) {
      if (!all(tocompare, [&](const typename PBCONTAINER2::value_type& pbcmp){ return fn(pbfilt, pbcmp); })) {
        selected += pbfilt;
      }
    }
    return selected;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void idiscardIfAll(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                            typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                        const typename PBCONTAINER2::value_type&)> fn) {
    tofilter = discardIfAll(tofilter, tocompare, fn);
  }


  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 selectIfAll(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                                  typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                              const typename PBCONTAINER2::value_type&)> fn) {
    PBCONTAINER1 selected;
    for (const auto& pbfilt : tofilter) {
      if (all(tocompare, [&](const typename PBCONTAINER2::value_type& pbcmp){ return fn(pbfilt, pbcmp); })) {
        selected += pbfilt;
      }
    }
    return selected;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void iselectIfAll(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare,
                           typename std::function<bool(const typename PBCONTAINER1::value_type&,
                                                       const typename PBCONTAINER2::value_type&)> fn) {
    tofilter = selectIfAll(tofilter, tocompare, fn);
  }

  /// @}


  /// @defgroup particlebaseutils_iso Isolation helpers
  /// @{

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void idiscardIfAnyDeltaRLess(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dR) {
    for (const typename PBCONTAINER2::value_type& pb : tocompare) {
      idiscard(tofilter, deltaRLess(pb, dR));
    }
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 discardIfAnyDeltaRLess(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dR) {
    PBCONTAINER1 tmp{tofilter};
    idiscardIfAnyDeltaRLess(tmp, tocompare, dR);
    return tmp;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void idiscardIfAnyDeltaPhiLess(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dphi) {
    for (const typename PBCONTAINER2::value_type& pb : tocompare) {
      idiscard(tofilter, deltaPhiLess(pb, dphi));
    }
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 discardIfAnyDeltaPhiLess(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dphi) {
    PBCONTAINER1 tmp{tofilter};
    idiscardIfAnyDeltaPhiLess(tmp, tocompare, dphi);
    return tmp;
  }



  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 selectIfAnyDeltaRLess(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dR) {
    PBCONTAINER1 selected;
    for (const typename PBCONTAINER1::value_type& f : tofilter) {
      if (any(tocompare, deltaRLess(f, dR))) selected.push_back(f);
    }
    return selected;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void iselectIfAnyDeltaRLess(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dR) {
    tofilter = selectIfAnyDeltaRLess(tofilter, tocompare, dR);
  }


  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline PBCONTAINER1 selectIfAnyDeltaPhiLess(const PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dphi) {
    PBCONTAINER1 selected;
    for (const typename PBCONTAINER1::value_type& f : tofilter) {
      if (any(tocompare, deltaPhiLess(f, dphi))) selected.push_back(f);
    }
    return selected;
  }

  template<typename PBCONTAINER1, typename PBCONTAINER2>
  inline void iselectIfAnyDeltaPhiLess(PBCONTAINER1& tofilter, const PBCONTAINER2& tocompare, double dphi) {
    tofilter = selectIfAnyDeltaPhiLess(tofilter, tocompare, dphi);
  }


  /// @todo Add 'all' variants

  /// @}


  /// @defgroup particlebaseutils_kin Unbound functions for kinematic properties
  ///
  /// @todo Mostly move to functions on FourMomentum
  /// @note In a sub-namespace (imported by default) for protection
  /// @{
  namespace Kin {

    /// Unbound function access to momentum
    inline FourMomentum mom(const ParticleBase& p) { return p.mom(); }
    /// Unbound function access to momentum
    inline FourMomentum p4(const ParticleBase& p) { return p.mom(); }

    /// Unbound function access to p3
    inline Vector3 p3(const ParticleBase& p) { return p.p3(); }

    /// Unbound function access to pTvec
    inline Vector3 pTvec(const ParticleBase& p) { return p.pTvec(); }

    /// Unbound function access to p
    inline double p(const ParticleBase& p) { return p.p(); }

    /// Unbound function access to pT
    inline double pT(const ParticleBase& p) { return p.pT(); }

    /// Unbound function access to E
    inline double E(const ParticleBase& p) { return p.E(); }

    /// Unbound function access to ET
    inline double Et(const ParticleBase& p) { return p.Et(); }

    /// Unbound function access to eta
    inline double eta(const ParticleBase& p) { return p.eta(); }

    /// Unbound function access to abseta
    inline double abseta(const ParticleBase& p) { return p.abseta(); }

    /// Unbound function access to rapidity
    inline double rap(const ParticleBase& p) { return p.rap(); }

    /// Unbound function access to abs rapidity
    inline double absrap(const ParticleBase& p) { return p.absrap(); }

    /// Unbound function access to mass
    inline double mass(const ParticleBase& p) { return p.mass(); }


    // /// Unbound function access to pair pT
    // inline double pairPt(const ParticleBase& p1, const ParticleBase& p2) { return (p1.mom() + p2.mom()).pT(); }

    // /// Unbound function access to pair mass
    // inline double pairMass(const ParticleBase& p1, const ParticleBase& p2) { return (p1.mom() + p2.mom()).mass(); }


    /// @brief Get the mass of a ParticleBase and a P4
    inline double mass(const ParticleBase& p, const FourMomentum& p4) {
      return mass(p.mom(), p4);
    }

    /// @brief Get the mass of a ParticleBase and a P4
    inline double mass(const FourMomentum& p4, const ParticleBase& p) {
      return mass(p4, p.mom());
    }

    /// @brief Get the mass of a pair of ParticleBase (as separate args)
    inline double mass(const ParticleBase& p1, const ParticleBase& p2) {
      return mass(p1.mom(), p2.mom());
    }

    /// @brief Get the mass^2 of a ParticleBase and a P4
    inline double mass2(const ParticleBase& p, const P4& p4) {
      return mass2(p.mom(), p4);
    }

    /// @brief Get the mass^2 of a ParticleBase and a P4
    inline double mass2(const P4& p4, const ParticleBase& p) {
      return mass2(p4, p.mom());
    }

    /// @brief Get the mass^2 of a pair of ParticleBase (as separate args)
    inline double mass2(const ParticleBase& p1, const ParticleBase& p2) {
      return mass2(p1.mom(), p2.mom());
    }

    /// @brief Get the transverse mass of a ParticleBase and a P4
    ///
    /// @note This ignores the particle mass and just computes mT from the 3-vectors
    /// @todo Fix!!
    inline double mT(const ParticleBase& p, const P4& p4) {
      return mT(p.mom(), p4);
    }

    /// @brief Get the transverse mass of a ParticleBase and a P4
    ///
    /// @note This ignores the particle mass and just computes mT from the 3-vectors
    /// @todo Fix!!
    inline double mT(const P4& p4, const ParticleBase& p) {
      return mT(p4, p.mom());
    }

    /// @brief Get the transverse mass of a pair of ParticleBase (as separate args)
    ///
    /// @note This ignores the particle mass and just computes mT from the 3-vectors
    /// @todo Fix!!
    inline double mT(const ParticleBase& p1, const ParticleBase& p2) {
      return mT(p1.mom(), p2.mom());
    }

    /// @brief Get the transverse momentum of a ParticleBase and a P4
    inline double pT(const ParticleBase& p, const P4& p4) {
      return pT(p.mom(), p4);
    }

    /// @brief Get the transverse momentum of a ParticleBase and a P4
    inline double pT(const P4& p4, const ParticleBase& p) {
      return pT(p4, p.mom());
    }

    /// @brief Get the transverse momentum of a pair of ParticleBase (as separate args)
    inline double pT(const ParticleBase& p1, const ParticleBase& p2) {
      return pT(p1.mom(), p2.mom());
    }

  }

  // Import Kin namespace into Rivet
  using namespace Kin;

  /// @}


  /// @defgroup pbcontcombutils ParticleBase container-combinatorics utils
  /// @{

  /// @brief Return the index from a vector which best matches mass(c[i]) to the target value
  ///
  /// A specialisation of closestMatchIndex from Utils.hh, with the
  /// function bound to Kin::mass as a common use-case.
  // auto closestMassIndex = std::bind(closestMatchIndex, std::placeholders::_1, Kin::mass, std::placeholders::_2, std::placeholders::_3);
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline int closestMassIndex(CONTAINER&& c, double mtarget, double mmin=-DBL_MAX, double mmax=DBL_MAX) {
    return closestMatchIndex(std::forward<CONTAINER>(c), Kin::mass, mtarget, mmin, mmax);
  }

  /// @brief Return the indices from two vectors which best match mass(c1[i], c2[j]) to the target value
  ///
  /// A specialisation of closestMatchIndex from Utils.hh, with the
  /// function bound to Kin::mass as a common use-case.
  template <typename CONTAINER1, typename CONTAINER2, typename = isCIterable<CONTAINER1, CONTAINER2>>
  inline pair<int,int> closestMassIndices(CONTAINER1&& c1, CONTAINER2&& c2,
        				   double mtarget, double mmin=-DBL_MAX, double mmax=DBL_MAX) {
    return closestMatchIndices(std::forward<CONTAINER1>(c1),
                               std::forward<CONTAINER2>(c2), Kin::mass, mtarget, mmin, mmax);
  }

  /// @brief Return the index from a vector which best matches mass(c[i], x) to the target value
  ///
  /// A specialisation of closestMatchIndex from Utils.hh, with the
  /// function bound to Kin::mass as a common use-case.
  template <typename CONTAINER, typename T, typename = isCIterable<CONTAINER>>
  inline int closestMassIndex(CONTAINER&& c, const T& x,
                              double mtarget, double mmin=-DBL_MAX, double mmax=DBL_MAX) {
    using FN = double(const ParticleBase&,const T&);
    return closestMatchIndex<CONTAINER, T, FN>(std::forward<CONTAINER>(c), x, Kin::mass, mtarget, mmin, mmax);
  }


  /// @brief Return the index from a vector which best matches mass(x, c[j]) to the target value
  ///
  /// A specialisation of closestMatchIndex from Utils.hh, with the
  /// function bound to Kin::mass as a common use-case.
  template <typename CONTAINER, typename T, typename = isCIterable<CONTAINER>>
  inline int closestMassIndex(T&& x, CONTAINER&& c,
                              double mtarget, double mmin=-DBL_MAX, double mmax=DBL_MAX) {
    return closestMatchIndex(std::forward<T>(x), std::forward<CONTAINER>(c), Kin::mass, mtarget, mmin, mmax);
  }

  /// @}


  /// @defgroup pbcontphysutils ParticleBase-container physics utils
  /// @{

  /// @brief Return the pT sum of a collection of particle-like objects
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline double sumPt(CONTAINER&& c) {
    return sum(c, Kin::pT, 0.0);
  }

  /// @brief Return the four-momentum sum of a collection of particle-like objects
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline FourMomentum sumP4(CONTAINER&& c) {
    return sum(c, Kin::p4, FourMomentum());
  }

  /// @brief Return the four-momentum sum of a collection of particle-like objects
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline Vector3 sumP3(CONTAINER&& c) {
    return sum(c, Kin::p3, Vector3());
  }

  /// @brief Return the four-momentum sum of a collection of particle-like objects
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline Vector3 sumPtVec(CONTAINER&& c) {
    return sum(c, Kin::pTvec, Vector3());
  }


  /// @brief Calculate the minimum phi separation in a collection of particle-like objects
  ///
  /// @note Returns in the range 0..pi
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline double deltaPhiMin(CONTAINER&& c) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c.size(); ++i) {
      for (size_t j = i; j < c.size(); ++j) {
        if (i == j) continue;
        const double dphi = fabs(deltaPhi(c[i], c[j]));
        if (dphi < rtn) rtn = dphi;
      }
    }
    return rtn;
  }

  /// @brief Calculate the minimum phi separation between two collections of particle-like objects
  ///
  /// @note Returns in the range 0..pi
  ///
  /// @warning Assumes the two vectors are mutually exclusive, so does not protect against self-comparisons.
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER1, typename CONTAINER2, typename = isCIterable<CONTAINER1, CONTAINER2>>
  inline double deltaPhiMin(CONTAINER1&& c1, CONTAINER2&& c2) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c1.size(); ++i) {
      for (size_t j = 0; j < c2.size(); ++j) {
        const double dphi = deltaPhi(c1[i], c2[j]);
        if (dphi < rtn) rtn = dphi;
      }
    }
    return rtn;
  }



  /// @brief Calculate the minimum pseudorapidity separation in a collection of particle-like objects
  ///
  /// @note Returns positive values
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline double deltaEtaMin(CONTAINER&& c) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c.size(); ++i) {
      for (size_t j = i; j < c.size(); ++j) {
        if (i == j) continue;
        const double deta = fabs(deltaEta(c[i], c[j]));
        if (deta < rtn) rtn = deta;
      }
    }
    return rtn;
  }

  /// @brief Calculate the minimum pseudorapidity separation between two collections of particle-like objects
  ///
  /// @note Returns positive value
  ///
  /// @warning Assumes the two vectors are mutually exclusive, so does not protect against self-comparisons.
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER1, typename CONTAINER2, typename = isCIterable<CONTAINER1, CONTAINER2>>
  inline double deltaEtaMin(CONTAINER1&& c1, CONTAINER2&& c2) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c1.size(); ++i) {
      for (size_t j = 0; j < c2.size(); ++j) {
        const double deta = deltaEta(c1[i], c2[j]);
        if (deta < rtn) rtn = deta;
      }
    }
    return rtn;
  }


  /// @brief Calculate the minimum rapidity separation in a collection of particle-like objects
  ///
  /// @note Returns positive values
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline double deltaRapMin(CONTAINER&& c) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c.size(); ++i) {
      for (size_t j = i; j < c.size(); ++j) {
        if (i == j) continue;
        const double dy = fabs(deltaRap(c[i], c[j]));
        if (dy < rtn) rtn = dy;
      }
    }
    return rtn;
  }

  /// @brief Calculate the minimum pseudorapidity separation between two collections of particle-like objects
  ///
  /// @note Returns positive values
  ///
  /// @warning Assumes the two vectors are mutually exclusive, so does not protect against self-comparisons.
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER1, typename CONTAINER2, typename = isCIterable<CONTAINER1, CONTAINER2>>
  inline double deltaRapMin(CONTAINER1&& c1, CONTAINER2&& c2) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c1.size(); ++i) {
      for (size_t j = 0; j < c2.size(); ++j) {
        const double dy = deltaRap(c1[i], c2[j]);
        if (dy < rtn) rtn = dy;
      }
    }
    return rtn;
  }



  /// @brief Calculate the minimum rapidity separation in a collection of particle-like objects
  ///
  /// @note Returns positive values
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER, typename = isCIterable<CONTAINER>>
  inline double deltaRMin(CONTAINER&& c) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c.size(); ++i) {
      for (size_t j = i; j < c.size(); ++j) {
        if (i == j) continue;
        const double dr = fabs(deltaR(c[i], c[j]));
        if (dr < rtn) rtn = dr;
      }
    }
    return rtn;
  }

  /// @brief Calculate the minimum pseudorapidity separation between two collections of particle-like objects
  ///
  /// @note Returns positive values
  ///
  /// @warning Assumes the two vectors are mutually exclusive, so does not protect against self-comparisons.
  ///
  /// @todo Generalise to min of any 2-ParticleBase function?
  template <typename CONTAINER1, typename CONTAINER2, typename = isCIterable<CONTAINER1, CONTAINER2>>
  inline double deltaRMin(CONTAINER1&& c1, CONTAINER2&& c2) {
    double rtn = DBL_MAX;
    for (size_t i = 0; i < c1.size(); ++i) {
      for (size_t j = 0; j < c2.size(); ++j) {
        const double dr = deltaR(c1[i], c2[j]);
        if (dr < rtn) rtn = dr;
      }
    }
    return rtn;
  }

  /// @}


  /// @}

}

#endif
