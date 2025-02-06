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
    HasBTag(const Cut& c=Cuts::open(), double dR=-1) : cut(c), deltaR(dR) {}
    // HasBTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.bTagged(cut, deltaR); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    const double deltaR;
  };
  using hasBTag = HasBTag;

  /// C-tagging functor, with a tag selection cut as the stored state
  struct HasCTag : BoolJetFunctor {
    HasCTag(const Cut& c=Cuts::open(), double dR=-1) : cut(c), deltaR(dR) {}
    // HasCTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.cTagged(cut, deltaR); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    const double deltaR;
  };
  using hasCTag = HasCTag;

  /// Tau-tagging functor, with a tag selection cut as the stored state
  struct HasTauTag : BoolJetFunctor {
    HasTauTag(const Cut& c=Cuts::open(), double dR=-1) : cut(c), deltaR(dR) {}
    // HasTauTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.tauTagged(cut, deltaR); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    const double deltaR;
  };
  using hasTauTag = HasTauTag;

  /// Anti-B/C-tagging functor, with a tag selection cut as the stored state
  struct HasNoTag : BoolJetFunctor {
    HasNoTag(const Cut& c=Cuts::open(), double dR=-1, bool quarktagsonly=false) : cut(c), deltaR(dR), qtagsonly(quarktagsonly) {}
    HasNoTag(const Cut& c=Cuts::open(), bool quarktagsonly=false) : HasNoTag(c, -1, quarktagsonly) {}
    // HasNoTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return !j.bTagged(cut, deltaR) && !j.cTagged(cut, deltaR) && (qtagsonly || !j.tauTagged(cut, deltaR)); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
    const double deltaR;
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


  /// @defgroup jetutils_trim Trimming operations for Rivet::Jets
  /// @{

  /// @brief Trim subjets with insufficient pT fraction
  ///
  /// Removes subjet constituents of pT less than @param frac*(jetpt) in-place.
  ///
  /// @note Mainly useful for reclustered jets.
  ///
  /// @warning Keeps all tags(), which might be incorrect, depending on use case.
  /// TODO: Fix this behaviour if necessary.
  Jet& itrimJetsFrac(Jet &jet, const double frac);

  /// @brief Trim subjets with insufficient pT fraction @param frac
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  ///
  /// @note Works on (almost) arbritrary containers of Jet
  ///
  /// TODO: Does not work on set<Jet> because of const restrictions of set. Fix if necessary.
  template <
    typename CONTAINER,
    typename = std::enable_if_t<
      is_citerable_v<CONTAINER>,
      Jet
    >
  >
  CONTAINER& itrimJetsFrac(CONTAINER &jets, const double frac, const JetSorter &sortFunc=cmpMomByPt){
    for ( Jet &jet : jets ) itrimJetsFrac(jet, frac);
    isortBy(jets, sortFunc);
    return jets;
  }

  template <typename T, typename U>
  std::map<T, U>& itrimJetsFrac(std::map<T, U> &jetMap, const double frac, const JetSorter &sortFunc=cmpMomByPt){
    for ( auto &item : jetMap ) itrimJetsFrac(item.second, frac, sortFunc);
    return jetMap;
  }

  /// @brief Trim subjets with insufficient pT fraction
  ///
  /// Removes subjet constituents of pT less than @param frac*(jetpt) out-of-place.
  ///
  /// @note Mainly useful for reclustered jets.
  ///
  /// @warning Keeps all tags(), which might be incorrect, depending on use case.
  /// TODO: Fix this behaviour if necessary.
  inline Jet trimJetsFrac(const Jet &jet, double frac){
    Jet rtn = jet;
    return itrimJetsFrac(rtn, frac);
  }

  /// @brief Trim subjets with insufficient pT fraction
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  ///
  /// @note Works on (almost) arbritrary containers of Jet
  ///
  /// TODO: Does not work on set<Jet> because of const restrictions of set. Fix if necessary.
  template <
    typename... Args, typename CONTAINER,
    typename = std::enable_if_t<
      is_citerable_v<CONTAINER>,
      Jet
    >
  >
  CONTAINER trimJetsFrac(const CONTAINER &jets, Args&&... args){
    CONTAINER rtn = jets;
    return itrimJetsFrac(rtn, std::forward<Args>(args)...);
  }

  template <typename T, typename U, typename... Args>
  std::map<T, U> trimJetsFrac(const std::map<T, U> &jetMap, Args&&... args){
    std::map<T, U> rtn = jetMap;
    return itrimJetsFrac(rtn, std::forward<Args>(args)...);
  }

  /// @}

  /// @defgroup jetutils_filt Operations to apply a FastJet::Filter to PseudoJet(s)
  /// @{

  /// @brief Apply given FastJet::Filter @param filter to PseudoJet @param pj in-place
  ///
  /// @note Can be any kind of filtering/ grooming algorithm, e.g. trimming.
  PseudoJet& ifilterPseudoJets(PseudoJet &pj, const fastjet::Filter &filter);

  /// @brief Apply given FastJet::Filter @param filter to container of PseudoJet @param pjs in-place
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  template <
    typename CONTAINER,
    typename = std::enable_if_t<
      is_citerable_v<CONTAINER>,
      PseudoJet
    >
  >
  CONTAINER& ifilterPseudoJets(CONTAINER &pjs, const fastjet::Filter &filter, const JetSorter &sortFunc=cmpMomByPt){
    for ( PseudoJet& pj : pjs ) ifilterPseudoJets(pj, filter);
    isortBy(pjs, sortFunc);
    return pjs;
  }

  /// @brief Apply given FastJet::Filter @param filter to map @param pjMap in-place
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  template <typename T, typename U, typename... Args>
  std::map<T, U>& ifilterPseudoJets(std::map<T, U> &pjMap, Args&&... args){
    for ( auto &item : pjMap ) ifilterPseudoJets(item.second, std::forward<Args>(args)...);
    return pjMap;
  }

  /// @brief Apply given FastJet::Filter @param filter to PseudoJet @param pj out-of-place
  inline PseudoJet filterPseudoJets(const PseudoJet &pj, const fastjet::Filter &filter){
    PseudoJet rtn = pj;
    return ifilterPseudoJets(rtn, filter);
  }

  /// @brief Apply given FastJet::Filter @param filter to container of PseudoJet @param pjs out-of-place
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  template <
    typename... Args, typename CONTAINER,
    typename = std::enable_if_t<
      is_citerable_v<CONTAINER>,
      PseudoJet
    >
  >
  CONTAINER filterPseudoJets(const CONTAINER &pjs, Args&&... args){
    CONTAINER rtn = pjs;
    return ifilterPseudoJets(rtn, std::forward<Args>(args)...);
  }

  /// @brief Apply given FastJet::Filter @param filter to map @param pjMap out-of-place
  ///
  /// @note Sorts by optional @param sortFunc, default is sorting by (descending) pT
  template <typename T, typename U, typename... Args>
  std::map<T, U> filterPseudoJets(const std::map<T, U> &pjMap, Args&&... args){
    std::map<T, U> rtn = pjMap;
    return ifilterPseudoJets(rtn, std::forward<Args>(args)...);
  }
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
