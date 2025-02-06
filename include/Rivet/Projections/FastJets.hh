// -*- C++ -*-
#ifndef RIVET_FastJets_HH
#define RIVET_FastJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/RivetFastJet.hh"

#include "fastjet/SISConePlugin.hh"
#include "fastjet/ATLASConePlugin.hh"
#include "fastjet/CMSIterativeConePlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"
#include "fastjet/TrackJetPlugin.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/VariableRPlugin.hh"

#include "Rivet/Projections/PxConePlugin.hh"
#include "Rivet/Tools/TypeTraits.hh"

namespace Rivet {

  typedef std::shared_ptr<fastjet::JetDefinition::Plugin> FJPluginPtr;

  /// Find jets using jet algorithms via the FastJet package
  class FastJets : public JetFinder {
  public:

    using JetFinder::operator =;


    /// @name Constructors etc.
    /// @{

    /// Constructor from a FastJet JetDefinition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : JetFinder(fsp, usemuons, useinvis), _jdef(jdef),
        _adef(adef), _cuts(Cuts::OPEN)
    {
      _initBase();
    }

    /// Constructor from a FastJet JetDefinition with optional Cut argument
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : JetFinder(fsp, usemuons, useinvis), _jdef(jdef),
        _adef(adef), _cuts(c)
    {
      _initBase();
    }

    /// JetDefinition-based constructor with reordered args for easier specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             fastjet::AreaDefinition* adef,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, jdef, usemuons, useinvis, adef)
    {    }

    /// JetDefinition-based constructor with Cut argument and reordered args for easier
    /// specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             const fastjet::JetDefinition& jdef,
             fastjet::AreaDefinition* adef,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, jdef, c, usemuons, useinvis, adef)
    {    }


    /// Native argument constructor, using FastJet alg/scheme enums.
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(type, rparameter, recom), usemuons, useinvis, adef)
    {    }

    /// Native argument constructor with Cut argument, using FastJet alg/scheme enums.
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(type, rparameter, recom), c, usemuons, useinvis, adef)
    {    }

    /// Native argument constructor with reordered args for easier specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             fastjet::AreaDefinition* adef,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, type, recom, rparameter, usemuons, useinvis, adef)
    {    }

    /// Native argument constructor with Cut argument and reordered args for easier specification of jet area definition
    ///
    /// @warning The AreaDefinition pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    FastJets(const FinalState& fsp,
             fastjet::JetAlgorithm type,
             fastjet::RecombinationScheme recom, double rparameter,
             fastjet::AreaDefinition* adef,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, type, recom, rparameter, c, usemuons, useinvis, adef)
    {    }

    /// @brief Explicitly pass in an externally-constructed plugin
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             FJPluginPtr plugin,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(plugin.get()), usemuons, useinvis, adef)
    {
      _plugin = plugin;
    }

    /// @brief Explicitly pass in an externally-constructed plugin
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             FJPluginPtr plugin,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : FastJets(fsp, fastjet::JetDefinition(plugin.get()), c, usemuons, useinvis, adef)
    {
      _plugin = plugin;
    }

    /// @brief Explicitly pass in an externally-constructed plugin, with reordered args for easier specification of jet area definition
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             FJPluginPtr plugin,
             fastjet::AreaDefinition* adef,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, plugin, usemuons, useinvis, adef)
    {    }

    /// @brief Explicitly pass in an externally-constructed plugin, with reordered args for easier specification of jet area definition
    ///
    /// @warning Provided plugin and area definition pointers must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             FJPluginPtr plugin,
             fastjet::AreaDefinition* adef,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE)
      : FastJets(fsp, plugin, c, usemuons, useinvis, adef)
    {    }

    /// @brief Convenience constructor using Rivet enums for most common jet algs (including some plugins).
    ///
    /// For the built-in algs, E-scheme recombination is used. For full control
    /// of FastJet built-in jet algs, use the constructors from native-args or a
    /// plugin pointer.
    ///
    /// @warning Provided area definition pointer must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             JetAlg alg, double rparameter=-1,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : JetFinder(fsp, usemuons, useinvis),
        _adef(adef), _cuts(Cuts::OPEN)
    {
      _initBase();
      _initJdef(alg, rparameter);
    }

    /// @brief Convenience constructor using Cut argument and Rivet enums for most common jet algs (including some plugins).
    ///
    /// For the built-in algs, E-scheme recombination is used. For full control
    /// of FastJet built-in jet algs, use the constructors from native-args or a
    /// plugin pointer.
    ///
    /// @warning Provided area definition pointer must be heap-allocated; Rivet will store/delete via a shared_ptr
    FastJets(const FinalState& fsp,
             JetAlg alg, double rparameter,
             const Cut& c,
             JetMuons usemuons=JetMuons::ALL,
             JetInvisibles useinvis=JetInvisibles::NONE,
             fastjet::AreaDefinition* adef=nullptr)
      : JetFinder(fsp, usemuons, useinvis), _adef(adef), _cuts(c)
    {
      _initBase();
      _initJdef(alg, rparameter);
    }



    /// Clone on the heap.
    RIVET_DEFAULT_PROJ_CLONE(FastJets);

    /// @}


    /// Import to avoid warnings about overload-hiding
    using Projection::operator =;


    /// @name Static helper functions for FastJet interaction, with tagging
    /// @{

    /// Make a FastJet JetDefinition according to enum JetAlg
    static fastjet::JetDefinition mkJetDef(JetAlg alg, double rparameter);

    /// Make a shared pointer to a FastJet plugin according to enum JetAlg
    /// Templated version for setting arbitrary parameters
    template<JetAlg JETALG, typename... Args>
    static FJPluginPtr mkPlugin(Args&&... args){
      return std::make_shared<mapJetAlg2Plugin_t<JETALG>>(std::forward<Args>(args)...);
    }

    /// Non-templated version only allowing to set rparameter
    static FJPluginPtr mkPlugin(JetAlg alg, double rparameter=-1);

    /// Make PseudoJets for input to a ClusterSequence, with user_index codes for constituent- and tag-particle linking
    static PseudoJets mkClusterInputs(const Particles& fsparticles, const Particles& tagparticles=Particles());
    /// Make a Rivet Jet from a PseudoJet holding a user_index code for lookup of Rivet fsparticle or tagparticle links
    static Jet mkJet(const PseudoJet& pj, const Particles& fsparticles, const Particles& tagparticles=Particles());
    /// Convert a whole list of PseudoJets to a list of Jets, with mkJet-style unpacking
    static Jets mkJets(const PseudoJets& pjs, const Particles& fsparticles=Particles(), const Particles& tagparticles=Particles());

    /// @}


    /// @defgroup jetutils_recluster QoL operations for FastJet reclustering
    /// @{

    typedef std::pair<PseudoJets, Particles> PJetsParts;

    /// @brief Recluster Rivet::Jets @param jetsIn according to fastjet::JetDefinition @param jDef
    ///
    /// @return jet collection with FastJet jet algorithm applied
    ///
    /// @note forwards the particles and tags of the original jets to the reclustered jets,
    /// i.e. R=0.8 jets reclustered from R=0.4 jets will still have the particles/tags of the R=0.4 jets as constituents,
    /// not the R=0.4 jets themselves as constituents.
    /// This means that the number of constituents (omitting the filter) is preserved
    ///
    /// @note Returned jets are sorted by (descending) pT
    template <
      typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, const fastjet::JetDefinition &jDef){
      PJetsParts reclusteredConsts = reclusterJetsParts(jetsIn, jDef);
      return mkTaggedJets(jetsIn, reclusteredConsts);
    }

    /// @brief Recluster Rivet::Jets @param jetsIn according to fastjet::JetDefinition @param jDef
    ///
    /// @return jet collection with FastJet jet algorithm applied
    ///
    /// @param filter a FastJet::Filter to be applied before re-assigning the tags
    /// Needed as separate function because FastJet doesn't grant us any easy initialisation check for filters
    ///
    /// @note Returned jets are sorted by (descending) pT
    template <
      typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, const fastjet::JetDefinition &jDef, const fastjet::Filter &filter){
      PJetsParts reclusteredConsts = reclusterJetsParts(jetsIn, jDef);
      ifilterPseudoJets(reclusteredConsts.first, filter);
      return mkTaggedJets(jetsIn, reclusteredConsts);
    }

    /// @brief Apply the @param jetAlg (given as shared pointer to JetDefinition::Plugin) to the given @param jetsIn jet collection.
    ///
    /// @return jet collection with FastJet jet algorithm applied
    template <
      typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, const FJPluginPtr &jetAlg){
      const fastjet::JetDefinition jDef(jetAlg.get());
      return reclusterJets(jetsIn, jDef);
    }

    /// @brief Apply the @param jetAlg (given as shared pointer to JetDefinition::Plugin) to the given @param jetsIn jet collection.
    ///
    /// @return jet collection with FastJet jet algorithm applied
    ///
    /// @param filter a FastJet::Filter to be applied before re-assigning the tags
    /// Needed as separate function because FastJet doesn't grant us any easy initialisation check for filters
    template <
      typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, const FJPluginPtr &jetAlg, const fastjet::Filter &filter){
      const fastjet::JetDefinition jDef(jetAlg.get());
      return reclusterJets(jetsIn, jDef, filter);
    }

    /// @brief Apply the @param JETALG (given as an template-enum JetAlg and arbitrary @param args...) to the given @param jetsIn jet collection.
    ///
    /// @return jet collection with FastJet jet algorithm applied
    template <
      JetAlg JETALG, typename... Args, typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, Args&&... args){
      if constexpr (JETALG<JetAlg::SISCONE){ // JetDefinition
        const fastjet::JetDefinition jDef = mkJetDef(JETALG, std::forward<Args>(args)...);
        return reclusterJets(jetsIn, jDef);
      } else { // Plugin
        const FJPluginPtr plugin = mkPlugin<JETALG>(std::forward<Args>(args)...);
        return reclusterJets(jetsIn, plugin);
      }
      throw std::invalid_argument( "Unknown jet algorithm: "+to_string(int(JETALG)) );
    }

    /// @brief Apply the @param JETALG (given as an template-enum JetAlg and arbitrary @param args...) to the given @param jetsIn jet collection.
    ///
    /// @return jet collection with FastJet jet algorithm applied
    ///
    /// @param filter a FastJet::Filter to be applied before re-assigning the tags
    /// Needed as separate function because FastJet doesn't grant us any easy initialisation check for filters
    template <
      JetAlg JETALG, typename... Args, typename CONTAINER,
      typename = std::enable_if_t<
        is_citerable_v<CONTAINER>,
        Jet
      >
    >
    static CONTAINER reclusterJets(const CONTAINER &jetsIn, const fastjet::Filter &filter, Args&&... args){
      if constexpr (JETALG<JetAlg::SISCONE){ // JetDefinition
        const fastjet::JetDefinition jDef = mkJetDef(JETALG, std::forward<Args>(args)...);
        return reclusterJets(jetsIn, jDef, filter);
      } else { // Plugin
        const FJPluginPtr plugin = mkPlugin<JETALG>(std::forward<Args>(args)...);
        return reclusterJets(jetsIn, plugin, filter);
      }
      throw std::invalid_argument( "Unknown jet algorithm: "+to_string(int(JETALG)) );
    }

    /// @brief Apply the @param args... to the given @param jetsMap map of jets.
    ///
    /// @return map of jets with FastJet jet algorithm applied
    template <typename T, typename U, typename... Args>
    static std::map<T, U> reclusterJets(const std::map<T, U> &jetsMap, Args&&... args){
      std::map<T, U> rtn;
      for ( auto const &[key, jetsIn] : jetsMap ) rtn[key] = reclusterJets(jetsIn, std::forward<Args>(args)...);
      return rtn;
    }

    /// @brief Apply the @param JETALG (given as an template-enum JetAlg and arbitrary @param args...) to the given @param jetsMap map of jets.
    ///
    /// @return map of jets with FastJet jet algorithm applied
    template <JetAlg JETALG, typename T, typename U, typename... Args>
    static std::map<T, U> reclusterJets(const std::map<T, U> &jetsMap, Args&&... args){
      std::map<T, U> rtn;
      for ( auto const &[key, jetsIn] : jetsMap ) rtn[key] = reclusterJets<JETALG>(jetsIn, std::forward<Args>(args)...);
      return rtn;
    }

    /// @}


    /// Reset the projection. Jet def, etc. are unchanged.
    void reset();


    /// @name Jet-area calculations
    /// @{

    /// @brief Use provided jet area definition
    ///
    /// @warning The provided pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    void useJetArea(fastjet::AreaDefinition* adef) {
      _adef.reset(adef);
    }

    /// Don't calculate a jet area
    void clearJetArea() {
      _adef.reset();
    }

    /// @}


    /// @name Jet grooming
    /// @{

    /// @brief Add a grooming transformer (base class of fastjet::Filter, etc.)
    ///
    /// @warning The provided pointer must be heap-allocated: it will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    void addTrf(fastjet::Transformer* trf) {
      _trfs.push_back(shared_ptr<fastjet::Transformer>(trf));
    }

    /// @brief Add a list of grooming transformers
    ///
    /// @warning The provided pointers must be heap-allocated: they will be stored/deleted via a shared_ptr.
    /// @note Provide an adef null pointer to re-disable jet area calculation
    template<typename TRFS, typename TRF=typename TRFS::value_type>
    typename std::enable_if<Derefable<TRF>::value, void>::type
    addTrfs(const TRFS& trfs) {
      for (auto& trf : trfs) addTrf(trf);
    }

    /// @brief Add a grooming filter
    ///
    /// @note Technically just a more physics-oriented wrapper of addTrf()
    void addFilter(fastjet::Filter* filter) {
      addTrf(filter);
    }

    /// @brief Add a list of grooming filters
    ///
    /// @note Technically just a more physics-oriented wrapper of addTrfs()
    template<typename FILTERS, typename FILTER=typename FILTERS::value_type>
    typename std::enable_if<Derefable<FILTER>::value, void>::type
    addFilters(const FILTERS& filters) {
      addTrfs(filters);
    }

    /// Don't apply any jet transformers
    void clearTrfs() {
      _trfs.clear();
    }

    /// @}


    /// @name Access to the jets
    /// @{

    /// Get the jets (unordered) with pT > ptmin.
    Jets _jets() const;

    /// Get the pseudo jets (unordered).
    /// @deprecated Use pseudojets
    PseudoJets pseudojets(double ptmin=0.0) const;

    /// Get the pseudo jets, ordered by \f$ p_T \f$.
    PseudoJets pseudojetsByPt(double ptmin=0.0) const {
      return sorted_by_pt(pseudojets(ptmin));
    }

    /// Get the pseudo jets, ordered by \f$ E \f$.
    PseudoJets pseudojetsByE(double ptmin=0.0) const {
      return sorted_by_E(pseudojets(ptmin));
    }

    /// Get the pseudo jets, ordered by rapidity.
    PseudoJets pseudojetsByRapidity(double ptmin=0.0) const {
      return sorted_by_rapidity(pseudojets(ptmin));
    }

    /// @}


    /// @name Access to the FastJet clustering objects such as jet def, area def, and cluster
    /// @{

    /// Return the cluster sequence.
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::ClusterSequence> clusterSeq() const {
      return _cseq;
    }

    /// Return the area-enabled cluster sequence (if an area defn exists, otherwise returns a null ptr).
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::ClusterSequenceArea> clusterSeqArea() const {
      return areaDef() ? dynamic_pointer_cast<fastjet::ClusterSequenceArea>(_cseq) : nullptr;
    }

    /// Return the jet definition.
    const fastjet::JetDefinition& jetDef() const {
      return _jdef;
    }

    /// @brief Return the area definition.
    ///
    /// @warning May be null!
    /// @todo Care needed re. const shared_ptr<T> vs. shared_ptr<const T>
    const shared_ptr<fastjet::AreaDefinition> areaDef() const {
      return _adef;
    }

    /// @}


  protected:

    /// Shared utility functions to implement constructor behaviour
    void _initBase();


    void _initJdef(JetAlg alg, double rparameter){
      MSG_DEBUG("JetAlg = " << static_cast<int>(alg));
      MSG_DEBUG("R parameter = " << rparameter);
      if ( alg < JetAlg::SISCONE ){ // fastjet::JetDefinitions
        _jdef = mkJetDef(alg, rparameter);
      } else { // fastjet::JetDefinition::Plugins
        _plugin = mkPlugin(alg, rparameter);
        _jdef = fastjet::JetDefinition(_plugin.get());
      }
    }

    /// Aux function for reclustering
    static Jets mkTaggedJets(const Jets &jetsIn, const PJetsParts &pJetsParts);

    /// Aux function for reclustering
    static PJetsParts reclusterJetsParts(const Jets &jetsIn, const fastjet::JetDefinition &jDef);

  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  public:

    /// Do the calculation locally (no caching).
    void calc(const Particles& fsparticles, const Particles& tagparticles=Particles());


  protected:

    /// Jet definition
    fastjet::JetDefinition _jdef;

    /// Pointer to user-handled area definition
    std::shared_ptr<fastjet::AreaDefinition> _adef;

    /// Cluster sequence
    std::shared_ptr<fastjet::ClusterSequence> _cseq;

    /// The kinematic cuts
    Cut _cuts;

    /// FastJet external plugin
    FJPluginPtr _plugin;

    /// List of jet groomers to be applied
    std::vector< std::shared_ptr<fastjet::Transformer> > _trfs;

    /// Map of Rivet::JetAlg to targeted fastjet::JetAlgorithm and fastjet::RecombinationScheme
    static const std::map<JetAlg, std::pair<fastjet::JetAlgorithm, fastjet::RecombinationScheme>> jetAlgMap;

    /// "Map" of Rivet::JetAlg to targeted fastjet::JetDefinition::Plugin
    /// Use dummy X so we can specialise mapJetAlg2Plugin in header
    template<JetAlg, typename X> struct mapJetAlg2Plugin;

    /// populate "map"
    template<typename X> struct mapJetAlg2Plugin<JetAlg::SISCONE,     X>{ using type = fastjet::SISConePlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::PXCONE,      X>{ using type = Rivet::PxConePlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::CDFJETCLU,   X>{ using type = fastjet::CDFJetCluPlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::CDFMIDPOINT, X>{ using type = fastjet::CDFMidPointPlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::D0ILCONE,    X>{ using type = fastjet::D0RunIIConePlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::JADE,        X>{ using type = fastjet::JadePlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::TRACKJET,    X>{ using type = fastjet::TrackJetPlugin; };
    template<typename X> struct mapJetAlg2Plugin<JetAlg::VARIABLER,   X>{ using type = fastjet::contrib::VariableRPlugin; };

    /// Make usage of "map" a bit more convenient
    template<JetAlg JETALG>
    using mapJetAlg2Plugin_t = typename mapJetAlg2Plugin<JETALG, void>::type;

    /// Map of vectors of y scales. This is mutable so we can use caching/lazy evaluation.
    mutable std::map<int, vector<double> > _yscales;

    /// Particles used for constituent and tag lookup
    Particles _fsparticles, _tagparticles;

  };

}

#endif
