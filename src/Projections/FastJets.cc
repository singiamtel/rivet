// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  void FastJets::_initBase() {
    setName("FastJets");
    declare(HeavyHadrons(), "HFHadrons");
    declare(TauFinder(TauDecay::HADRONIC), "Taus");

    // Print/hide FJ banner
    std::cout.setstate(std::ios_base::badbit);
    fastjet::ClusterSequence::print_banner();
    std::cout.clear();
  }

  const std::map<JetAlg, std::pair<fastjet::JetAlgorithm, fastjet::RecombinationScheme>> FastJets::jetAlgMap = {
    {JetAlg::KT,        std::make_pair(fastjet::kt_algorithm,        fastjet::E_scheme)},
    {JetAlg::ANTIKT,    std::make_pair(fastjet::antikt_algorithm,    fastjet::E_scheme)},
    {JetAlg::CAM,       std::make_pair(fastjet::cambridge_algorithm, fastjet::E_scheme)},
    {JetAlg::DURHAM,    std::make_pair(fastjet::ee_kt_algorithm,     fastjet::E_scheme)},
    {JetAlg::GENKTEE,   std::make_pair(fastjet::ee_genkt_algorithm,  fastjet::external_scheme)},
    {JetAlg::KTET,      std::make_pair(fastjet::kt_algorithm,        fastjet::Et_scheme)},
    {JetAlg::ANTIKTET,  std::make_pair(fastjet::antikt_algorithm,    fastjet::Et_scheme)}
  };

  fastjet::JetDefinition FastJets::mkJetDef(const JetAlg alg, const double rparameter){
    if (jetAlgMap.find(alg)==jetAlgMap.end()) {
      throw std::invalid_argument( to_string(int(alg)) + " is no known jet algorithm." );
    }
    const std::pair<fastjet::JetAlgorithm, fastjet::RecombinationScheme> p = jetAlgMap.at(alg);
    if (alg == JetAlg::DURHAM)  return fastjet::JetDefinition(p.first, p.second);
    return fastjet::JetDefinition(p.first, rparameter, p.second);
  }

  FJPluginPtr FastJets::mkPlugin(const JetAlg alg, const double rparameter){
    /// Unfortunately, we can't simply call the templated version of mkPlugin<> here because `alg` is not a constant expression
    if (alg == JetAlg::SISCONE) {
      const double OVERLAP_THRESHOLD = 0.75;
      return make_shared<fastjet::SISConePlugin>(rparameter, OVERLAP_THRESHOLD);
    }
    if (alg == JetAlg::PXCONE) {
      return make_shared<Rivet::PxConePlugin>(rparameter);
    }
    if (alg == JetAlg::CDFJETCLU) {
      const double OVERLAP_THRESHOLD = 0.75;
      const double SEED_THRESHOLD = 1.0;
      return make_shared<fastjet::CDFJetCluPlugin>(rparameter, OVERLAP_THRESHOLD, SEED_THRESHOLD);
    }
    if (alg == JetAlg::CDFMIDPOINT) {
      const double OVERLAP_THRESHOLD = 0.5;
      const double SEED_THRESHOLD = 1.0;
      return make_shared<fastjet::CDFMidPointPlugin>(rparameter, OVERLAP_THRESHOLD, SEED_THRESHOLD);
    }
    if (alg == JetAlg::D0ILCONE) {
      const double MIN_JET_ET = 6.0;
      return make_shared<fastjet::D0RunIIConePlugin>(rparameter, MIN_JET_ET);
    }
    if (alg == JetAlg::JADE) {
      return make_shared<fastjet::JadePlugin>();
    }
    if (alg == JetAlg::TRACKJET) {
      return make_shared<fastjet::TrackJetPlugin>(rparameter);
    }
    if (alg == JetAlg::VARIABLER) {
      const double RHO = 550*GeV;
      const double RMAX = 1.5;
      return make_shared<fastjet::contrib::VariableRPlugin>(RHO, rparameter, RMAX, fastjet::contrib::VariableRPlugin::AKTLIKE);
    }
    throw std::invalid_argument( "JetAlg " + to_string(int(alg)) + "cannot be converted to fastjet::JetDefinition::Plugin" );
  }

  CmpState FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    CmpState rtn =
      cmp(_useMuons, other._useMuons) ||
      cmp(_useInvisibles, other._useInvisibles) ||
      mkNamedPCmp(other, "FS") ||
      cmp(_jdef.jet_algorithm(), other._jdef.jet_algorithm()) ||
      cmp(_jdef.recombination_scheme(), other._jdef.recombination_scheme()) ||
      cmp(_jdef.plugin(), other._jdef.plugin()) ||
      cmp(_jdef.R(), other._jdef.R()) ||
      cmp(_cuts, other._cuts) || cmp(_adef, other._adef);
    if (rtn != CmpState::EQ) return rtn; //< shortcut transformer comparison if aleady different

    // Compare the transformers list
    if (_trfs.empty() && other._trfs.empty()) return CmpState::EQ;
    /// @todo Improve fastjet::Transformer to add a virtual operator==, and use all()
    // if (_trfs.size() != other._trfs.size()) return CmpState::NEQ;
    // for (size_t it = 0; it < _trfs.size(); ++it) {

    return CmpState::NEQ;
  }


  // STATIC
  PseudoJets FastJets::mkClusterInputs(const Particles& fsparticles, const Particles& tagparticles) {
    PseudoJets pjs;
    /// @todo Use FastJet3's UserInfo system to store Particle pointers directly?

    /// @todo Support DeltaR tagging in place of ghost association

    // Store 4 vector data about each particle into FastJet's PseudoJets
    for (size_t i = 0; i < fsparticles.size(); ++i) {
      fastjet::PseudoJet pj = fsparticles[i];
      pj.set_user_index(i+1);
      pjs.push_back(pj);
    }
    // And the same for ghost tagging particles (with negative user indices)
    for (size_t i = 0; i < tagparticles.size(); ++i) {
      fastjet::PseudoJet pj = tagparticles[i];
      pj *= 1e-20; ///< Ghostify the momentum
      pj.set_user_index(-i-1);
      pjs.push_back(pj);
    }

    return pjs;
  }


  // STATIC
  Jet FastJets::mkJet(const PseudoJet& pj, const Particles& fsparticles, const Particles& tagparticles) {
    const PseudoJets pjconstituents = pj.constituents();

    Particles constituents, tags;
    constituents.reserve(pjconstituents.size());

    /// @todo Support DeltaR tagging in place of ghost association

    for (const fastjet::PseudoJet& pjc : pjconstituents) {
      // Pure ghosts don't have corresponding particles
      if (pjc.has_area() && pjc.is_pure_ghost()) continue;
      // Default user index = 0 doesn't give valid particle lookup
      if (pjc.user_index() == 0) continue;
      // Split by index sign into constituent & tag lookup
      if (pjc.user_index() > 0) {
        // Find constituents if index > 0
        const size_t i = pjc.user_index() - 1;
        if (i >= fsparticles.size()) throw RangeError("FS particle lookup failed in jet construction");
        constituents.push_back(fsparticles.at(i));
      } else if (!tagparticles.empty()) {
        // Find tags if index < 0
        const size_t i = abs(pjc.user_index()) - 1;
        if (i >= tagparticles.size()) throw RangeError("Tag particle lookup failed in jet construction");
        tags.push_back(tagparticles.at(i));
      } else {
        // Treat as particle dummy
        constituents.push_back(Particle(0, Jet(pjc).mom()));
      }
    }

    return Jet(pj, constituents, tags);
  }


  // STATIC
  Jets FastJets::mkJets(const PseudoJets& pjs, const Particles& fsparticles, const Particles& tagparticles) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs) {
      rtn.push_back(FastJets::mkJet(pj, fsparticles, tagparticles));
    }
    return rtn;
  }

  FastJets::PJetsParts FastJets::reclusterJetsParts(const Jets &jetsIn, const fastjet::JetDefinition &jDef){
    // Collate jets as pseudoparticles for tagging
    Particles jetsParts;
    for(const Jet &j : jetsIn){
      // Add jet as pseudo-particle
      jetsParts.push_back(Particle(0, j.mom()));
    }
    PseudoJets pjsIn = FastJets::mkClusterInputs(jetsParts);

    // Recluster jets
    ClusterSequence *cSeq = new ClusterSequence(pjsIn, jDef);
    PseudoJets pjsReclustered = cSeq->inclusive_jets();
    cSeq->delete_self_when_unused();
    // Should delete pointer when all the Pseudojets associated with it have gone out of scope:
    // https://fastjet.fr/repo/doxygen-3.4.3/classfastjet_1_1ClusterSequence.html#a3a5f949d32342cc4bde9f3bcb33bbe47
    // Fingers crossed FastJet actually does that correctly and this doesn't cause a memory leak

    return PJetsParts(pjsReclustered, jetsParts);
  }

  Jets FastJets::mkTaggedJets(const Jets &jetsIn, const PJetsParts &pJetsParts){
    const PseudoJets &pjsReclustered = pJetsParts.first;

    // Convert to Rivet::Jets
    Jets jetsOut = FastJets::mkJets(pjsReclustered, pJetsParts.second, Particles());

    // Replace particles and tags of reclustered jets with those from original jets
    for ( size_t pos=0; pos<pjsReclustered.size(); ++pos ){
      const PseudoJet &pjLargeR = pjsReclustered[pos];
      Particles &particles = jetsOut[pos].particles();
      Particles &tags = jetsOut[pos].tags();

      // Reset particles and tags
      particles = Particles();
      tags = Particles();

      // Forward particles and tags from original jets
      for ( const PseudoJet &pjSmallR : pjLargeR.constituents() ){
        const Jet &jSmallR = jetsIn[pjSmallR.user_index()-1];

        particles.insert(particles.end(), jSmallR.particles().begin(), jSmallR.particles().end());
        tags.insert(tags.end(), jSmallR.tags().begin(), jSmallR.tags().end());
      }
    }

    // Sort in-place by pT, probably wanted in most cases
    isortByPt(jetsOut);

    return jetsOut;
  }

  void FastJets::project(const Event& e) {
    // Assemble final state particles
    const string fskey = (_useInvisibles == JetInvisibles::NONE) ? "VFS" : "FS";
    Particles fsparticles = apply<FinalState>(e, fskey).particles();
    // Remove prompt invisibles if needed (already done by VFS if using NO_INVISIBLES)
    if (_useInvisibles == JetInvisibles::DECAY) {
      idiscard(fsparticles, [](const Particle& p) { return !p.isVisible() && p.isPrompt(); });
    }
    // Remove prompt/all muons if needed
    if (_useMuons == JetMuons::DECAY) {
      idiscard(fsparticles, [](const Particle& p) { return isMuon(p) && p.isPrompt(); });
    } else if (_useMuons == JetMuons::NONE) {
      idiscard(fsparticles, isMuon);
    }

    // Tagging particles
    /// @todo Support DeltaR tagging in place of ghost association
    const Particles chadrons = apply<HeavyHadrons>(e, "HFHadrons").cHadrons();
    const Particles bhadrons = apply<HeavyHadrons>(e, "HFHadrons").bHadrons();
    const Particles taus = apply<FinalState>(e, "Taus").particles();

    // Run the calculation
    calc(fsparticles, chadrons+bhadrons+taus);
  }


  void FastJets::calc(const Particles& fsparticles, const Particles& tagparticles) {
    MSG_DEBUG("Finding jets from " << fsparticles.size() << " input particles + " << tagparticles.size() << " tagging particles");
    _fsparticles = fsparticles;
    _tagparticles = tagparticles;

    // Make pseudojets, with mapping info to Rivet FS and tag particles
    PseudoJets pjs = mkClusterInputs(_fsparticles, _tagparticles);

    // Run either basic or area-calculating cluster sequence as reqd.
    /// @todo How can we make sure these persist if we call calc multiple times?
    if (_adef) {
      _cseq.reset(new fastjet::ClusterSequenceArea(pjs, _jdef, *_adef));
    } else {
      _cseq.reset(new fastjet::ClusterSequence(pjs, _jdef));
    }

    MSG_DEBUG("ClusterSequence constructed; Njets_tot = "
              << _cseq->inclusive_jets().size() << ", Njets(pT > 10 GeV) = "
              << _cseq->inclusive_jets(10*GeV).size());
  }


  void FastJets::reset() {
    _yscales.clear();
    _fsparticles.clear();
    _tagparticles.clear();
    /// @todo _cseq = fastjet::ClusterSequence();
  }


  Jets FastJets::_jets() const {
    /// @todo Cache?
    return mkJets(pseudojets(), _fsparticles, _tagparticles);
  }


  PseudoJets FastJets::pseudojets(double ptmin) const {
    // Get the base set of pseudo-jets
    PseudoJets tmp = clusterSeq() ? clusterSeq()->inclusive_jets(ptmin) : PseudoJets();

    // Run the jet groomers on each jet
    PseudoJets rtn;
    rtn.reserve(ceil(0.5*tmp.size())); // anticipate half of them to pass on average
    for (PseudoJet& pj : tmp) {
      for (auto& t : _trfs) {
        pj = t->result(pj);
      }
      if (_cuts->accept(Jet(pj)))  rtn.push_back(std::move(pj));
    }

    return rtn;
  }


}
