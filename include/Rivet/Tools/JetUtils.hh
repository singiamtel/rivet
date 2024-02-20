#ifndef RIVET_JETUTILS_HH
#define RIVET_JETUTILS_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"

namespace Rivet {


  /// @defgroup jetutils Functions for Jets
  /// @{

  /// @defgroup jetutils_conv Converting between Jets, Particles and PseudoJets
  /// @{

  inline PseudoJets mkPseudoJets(const Particles& ps) {
    PseudoJets rtn; rtn.reserve(ps.size());
    for (const Particle& p : ps) rtn.push_back(p.pseudojet());
    return rtn;
  }

  inline PseudoJets mkPseudoJets(const Jets& js) {
    PseudoJets rtn; rtn.reserve(js.size());
    for (const Jet& j : js) rtn.push_back(j.pseudojet());
    return rtn;
  }

  inline Jets mkJets(const PseudoJets& pjs) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs) rtn.push_back(pj);
    return rtn;
  }

  /// @}


  /// @defgroup jetutils_j2bool Jet classifier -> bool functors
  /// @{

  /// std::function instantiation for functors taking a Jet and returning a bool
  using JetSelector = function<bool(const Jet&)>;
  /// std::function instantiation for functors taking two Jets and returning a bool
  using JetSorter = function<bool(const Jet&, const Jet&)>;


  /// Base type for Jet -> bool functors
  struct BoolJetFunctor {
    virtual bool operator()(const Jet& p) const = 0;
    virtual ~BoolJetFunctor() {}
  };


  /// Functor for and-combination of selector logic
  struct BoolJetAND : public BoolJetFunctor {
    BoolJetAND(const std::vector<JetSelector>& sels) : selectors(sels) {}
    BoolJetAND(const JetSelector& a, const JetSelector& b) : selectors({a,b}) {}
    BoolJetAND(const JetSelector& a, const JetSelector& b, const JetSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Jet& j) const {
      for (const JetSelector& sel : selectors) if (!sel(j)) return false;
      return true;
    }
    std::vector<JetSelector> selectors;
  };
  /// Operator syntactic sugar for AND construction
  inline BoolJetAND operator && (const JetSelector& a, const JetSelector& b) {
    return BoolJetAND(a, b);
  }


  /// Functor for or-combination of selector logic
  struct BoolJetOR : public BoolJetFunctor {
    BoolJetOR(const std::vector<JetSelector>& sels) : selectors(sels) {}
    BoolJetOR(const JetSelector& a, const JetSelector& b) : selectors({a,b}) {}
    BoolJetOR(const JetSelector& a, const JetSelector& b, const JetSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Jet& j) const {
      for (const JetSelector& sel : selectors) if (sel(j)) return true;
      return false;
    }
    std::vector<JetSelector> selectors;
  };
  /// Operator syntactic sugar for OR construction
  inline BoolJetOR operator || (const JetSelector& a, const JetSelector& b) {
    return BoolJetOR(a, b);
  }


  /// Functor for inverting selector logic
  struct BoolJetNOT : public BoolJetFunctor {
    BoolJetNOT(const JetSelector& sel) : selector(sel) {}
    bool operator()(const Jet& j) const { return !selector(j); }
    JetSelector selector;
  };
  /// Operator syntactic sugar for NOT construction
  inline BoolJetNOT operator ! (const JetSelector& a) {
    return BoolJetNOT(a);
  }



  /// B-tagging functor, with a tag selection cut as the stored state
  struct HasBTag : BoolJetFunctor {
    HasBTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasBTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.bTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasBTag = HasBTag;

  /// C-tagging functor, with a tag selection cut as the stored state
  struct HasCTag : BoolJetFunctor {
    HasCTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasCTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.cTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasCTag = HasCTag;

  /// Tau-tagging functor, with a tag selection cut as the stored state
  struct HasTauTag : BoolJetFunctor {
    HasTauTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasTauTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.tauTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasTauTag = HasTauTag;

  /// Anti-B/C-tagging functor, with a tag selection cut as the stored state
  struct HasNoTag : BoolJetFunctor {
    HasNoTag(const Cut& c=Cuts::open(), bool quarktagsonly=false) : cut(c), qtagsonly(quarktagsonly) {}
    // HasNoTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return !j.bTagged(cut) && !j.cTagged(cut) && (qtagsonly || !j.tauTagged(cut)); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    bool qtagsonly;
  };
  using hasNoTag = HasNoTag;

  /// @}


  /// @defgroup jetutils_filt Unbound functions for filtering jets
  /// @{

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  Jets& iselect(Jets& jets, const Cut& c);


  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets select(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return iselect(rtn, c);
  }


  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets select(const Jets& jets, const Cut& c, Jets& out) {
    out = select(jets, c);
    return out;
  }



  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  Jets& idiscard(Jets& jets, const Cut& c);


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets discard(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return idiscard(rtn, c);
  }


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets discard(const Jets& jets, const Cut& c, Jets& out) {
    out = discard(jets, c);
    return out;
  }

  /// @}


  /// @defgroup jetutils_trim Trimming operations not covered in fastjet
  /// @{

  /// Take in PseudoJets, Return Jets with subjet constituents of pt under
  /// frac*(jetpt) removed.
  /// Mainly useful for reclustered jets.
  Jets trimJetsFrac(const PseudoJets& jetsIn, const double frac=0.1);

  /// @}



  /// @defgroup jetutils_coll Operations on collections of Jet
  /// @note This can't be done on generic collections of ParticleBase -- thanks, C++ :-/
  /// @{
  namespace Kin {

    /// @todo This shouldn't be necessary, if the sum() function SFINAE picked up be ParticleBase versions...
    inline double pT(const Jet& j) {
      return j.pT();
    }

    inline double sumPt(const Jets& js) {
      return sum(js, Kin::pT, 0.0);
    }

    inline FourMomentum sumP4(const Jets& js) {
      return sum(js, Kin::p4, FourMomentum());
    }

    inline Vector3 sumP3(const Jets& js) {
      return sum(js, Kin::p3, Vector3());
    }

    /// @todo Min dPhi, min dR?
    /// @todo Isolation routines?

  }
  /// @}


  // Import Kin namespace into Rivet
  using namespace Kin;


  /// @}

}

#endif
