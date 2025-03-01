// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief Simplified analysis API with predefined physics objects
  class SimpleAnalysis : public Analysis {
  public:

    /// Physics-object ID classes
    enum class IDClass { LOOSE, MEDIUM, TIGHT };

    /// Constructor, passing arg to base class and supporting subclassing
    SimpleAnalysis(const std::string& name)
      : Analysis(name)
    {    }

    /// Destructor to support subclassing
    virtual ~SimpleAnalysis()
    {    }


    /// @name Analysis methods
    /// @{

    /// Initialise built-in projections before the run, after user init()
    ///
    /// @todo Only register default proj if the name hasn't already
    /// been used, to allow user overrides
    void postInit() {

      // Electron projections
      FinalState es(_acut && _ecut && Cuts::abspid == PID::ELECTRON);
      LeptonFinder es_prompt(_acut && _ecut && Cuts::abspid == PID::ELECTRON, 0.1);
      for (const auto& classkv : _idclasses()) {
        const IDClass& idc = classkv.first;
        const string p = "Electrons";
        const string aname = _pckey(true, idc, p);
        if (!hasProjection(aname))
          declare(SmearedParticles(es, _effs_electron[idc], _smear_electron), aname);
        const string pname = _pckey(false, idc, p);
        if (!hasProjection(pname))
          declare(SmearedParticles(es_prompt, _effs_electron[idc], _smear_electron), pname);
      }

      // Muon projections
      FinalState ms(_acut && _mcut && Cuts::abspid == PID::MUON);
      LeptonFinder ms_prompt(_acut && _mcut && Cuts::abspid == PID::MUON, 0.1);
      for (const auto& classkv : _idclasses()) {
        const IDClass& idc = classkv.first;
        const string p = "Muons";
        const string aname = _pckey(true, idc, p);
        if (!hasProjection(aname))
          declare(SmearedParticles(ms, _effs_muon[idc], _smear_muon), aname);
        const string pname = _pckey(false, idc, p);
        if (!hasProjection(pname))
          declare(SmearedParticles(ms_prompt, _effs_muon[idc], _smear_muon), pname);
      }

      // Photon projections
      FinalState ys(_acut && _phocut && Cuts::abspid == PID::PHOTON);
      DirectFinalState ys_prompt(_acut && _phocut && Cuts::abspid == PID::PHOTON);
      for (const auto& classkv : _idclasses()) {
        const IDClass& idc = classkv.first;
        const string p = "Photons";
        const string aname = _pckey(true, idc, p);
        if (!hasProjection(aname))
          declare(SmearedParticles(ys, _effs_photon[idc], _smear_photon), aname);
        const string pname = _pckey(false, idc, p);
        if (!hasProjection(pname))
          declare(SmearedParticles(ys_prompt, _effs_photon[idc], _smear_photon), pname);
      }

      // Tau projections
      TauFinder ts(TauDecay::HADRONIC, LeptonOrigin::ANY, _acut && _tcut);
      TauFinder ts_prompt(TauDecay::HADRONIC, LeptonOrigin::NODECAY, _acut && _tcut);
      for (const auto& classkv : _idclasses()) {
        const IDClass& idc = classkv.first;
        const string p = "Taus";
        const string aname = _pckey(true, idc, p);
        if (!hasProjection(aname))
          declare(SmearedParticles(ts, _effs_tau[idc], _smear_tau), aname);
        const string pname = _pckey(false, idc, p);
        if (!hasProjection(pname))
          declare(SmearedParticles(ts_prompt, _effs_tau[idc], _smear_tau), pname);
      }

      // Track projection
      const FinalState trkfs(_acut && _trkcut && Cuts::abscharge > 0);
      if (!hasProjection("Tracks"))
        declare(SmearedParticles(trkfs, _eff_trk, _smear_trk), "Tracks");

      // Jet projections
      const FinalState allfs(_acut);
      declare(allfs, "AllParticles");
      for (const auto& jditem : _jdefs) {
        const JetScheme& jd = jditem.second;
        /// @todo Allow multiple jet-smearing / tagging fns... how to default?
        FastJets jetfs(allfs, jd.alg, jd.R, jd.muons, jd.invis);
        /// @todo The override is more complicated here since users already control the names... throw an exception instead?
        if (!hasProjection("Jets"+jditem.first))
          declare(SmearedJets(jetfs, _smears_jet[jditem.first], _effs_btag[jditem.first]), "Jets"+jditem.first);
      }

      // Missing momentum projection
      if (!hasProjection("MET"))
        declare(SmearedMET(MissingMomentum(allfs), _smearps_met, _smear_met), "MET");
    }


    /// Reset per-event members, before the user analyze()
    void preAnalyze(const Event& event) {
      _electrons.clear();
      _muons.clear();
      _taus.clear();
      _photons.clear();
      _jets.clear();
      _bjets.clear();
    }

    /// @}



    /// @name Detector and reco config
    /// @{

    /// @name Detector smearing setup (call in init)
    /// @{

    /// Set the electron reco functions for a given ID class
    void setElectronReco(const ParticleEffFn& eff, const ParticleSmearFn& smear=PARTICLE_SMEAR_IDENTITY, IDClass idc=IDClass::MEDIUM) {
      _effs_electron[idc] = eff;
      _smear_electron = smear;
    }
    /// Set the muon reco functions for a given ID class
    void setMuonReco(const ParticleEffFn& eff, const ParticleSmearFn& smear=PARTICLE_SMEAR_IDENTITY, IDClass idc=IDClass::MEDIUM) {
      _effs_muon[idc] = eff;
      _smear_muon = smear;
    }
    /// Set the photon reco functions for a given ID class
    void setPhotonReco(const ParticleEffFn& eff, const ParticleSmearFn& smear=PARTICLE_SMEAR_IDENTITY, IDClass idc=IDClass::MEDIUM) {
      _effs_photon[idc] = eff;
      _smear_photon = smear;
    }
    /// Set the tau reco functions for a given ID class
    void setTauReco(const ParticleEffFn& eff, const ParticleSmearFn& smear=PARTICLE_SMEAR_IDENTITY, IDClass idc=IDClass::MEDIUM) {
      _effs_tau[idc] = eff;
      _smear_tau = smear;
    }
    /// Set the track reco functions for a given jet name
    void setTrkReco(const ParticleEffFn& eff, const ParticleSmearFn& smear=PARTICLE_SMEAR_IDENTITY) {
      _eff_trk = eff;
      _smear_trk = smear;
    }
    /// Set the jet reco functions for a given jet name
    void setJetReco(const JetSmearFn& smear, const JetEffFn& btageff, const std::string& jetname="DEFAULT") {
      _smears_jet[jetname] = smear;
      _effs_btag[jetname] = btageff;
    }
    /// Set the MET reco-params function
    void setMETReco(const METSmearParamsFn& smearps) {
      _smearps_met = smearps;
    }
    /// Set the MET reco functions
    void setMETReco(const METSmearParamsFn& smearps, const METSmearFn& smear) {
      _smearps_met = smearps;
      _smear_met = smear;
    }

    /// @}


    /// @name Baseline object config (call in init)
    ///
    /// @todo Allow specialising cuts to ID class
    /// @{
    void setAccCut(const Cut& cut) { _acut = cut; }
    void setElectronCut(const Cut& cut) { _ecut = cut; }
    void setMuonCut(const Cut& cut) { _mcut = cut; }
    void setPhotonCut(const Cut& cut) { _phocut = cut; }
    void setTauCut(const Cut& cut) { _tcut = cut; }
    void setTrkCut(const Cut& cut) { _trkcut = cut; }
    void setJetDef(const JetScheme& jdef, const std::string& jetname="DEFAULT") { _jdefs[jetname] = jdef; }
    void setJetCut(const Cut& cut, const std::string& jetname="DEFAULT") { _jcuts[jetname] = cut; }
    /// @}

    /// @}



    /// @name Obtaining physics objects
    /// @{

    /// Get the pT-ordered set of reconstructed prompt electrons
    const Particles& electrons(IDClass eclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      const std::string kpc = _pckey(include_nonprompt, eclass);
      if (_electrons.find(kpc) == _electrons.end()) { //< if not cached
        _electrons[kpc] = apply<ParticleFinder>(currentEvent(), "Electrons"+kpc).particlesByPt();
      }
      return _electrons.at(kpc);
    }
    /// Get the pT-ordered set of reconstructed prompt electrons, with extra cuts
    Particles electrons(const Cut& cut, IDClass eclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      return select(electrons(eclass, include_nonprompt), cut);
    }


    /// Get the pT-ordered set of reconstructed muons
    const Particles& muons(IDClass muclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      const std::string kpc = _pckey(include_nonprompt, muclass);
      if (_muons.find(kpc) == _muons.end()) { //< if not cached
        _muons[kpc] = apply<ParticleFinder>(currentEvent(), "Muons"+kpc).particlesByPt();
      }
      return _muons.at(kpc);
    }
    /// Get the pT-ordered set of reconstructed muons, with extra cuts
    Particles muons(const Cut& cut, IDClass muclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      return select(muons(muclass, include_nonprompt), cut);
    }


    /// Get the pT-ordered set of reconstructed prompt and isolated photons
    const Particles& photons(IDClass phoclass=IDClass::MEDIUM, bool include_nonprompt=false, double dRiso=0.4, double isofrac=0.1) const {
      const std::string kpc = _pckey(include_nonprompt, phoclass);
      if (_photons.find(kpc) == _photons.end()) { //< if not cached
        Particles ystmp = apply<ParticleFinder>(currentEvent(), "Photons"+kpc).particlesByPt();
        const Particles& clusters = apply<ParticleFinder>(currentEvent(), "AllParticles").particles();
        _photons[kpc].reserve(ystmp.size());
        for (const Particle& y : ystmp) {
          double sumpt = 0;
          for (const Particle& cl : select(clusters, deltaRLess(y, dRiso))) {
            if (deltaR(cl, y) == 0.01) continue; //< irresolvably close
            sumpt += cl.pT();
          }
          if (sumpt/y.pT() < isofrac) _photons[kpc].push_back(y);
        }
      }
      return _photons.at(kpc);
    }
    /// Get the pT-ordered set of reconstructed prompt and isolated photons, with extra cuts
    Particles photons(const Cut& cut, IDClass phoclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      return select(photons(phoclass, include_nonprompt), cut);
    }


    /// Get the pT-ordered set of reconstructed prompt (hadronic) taus
    const Particles& taus(IDClass tauclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      const std::string kpc = _pckey(include_nonprompt, tauclass);
      if (_taus.find(kpc) == _taus.end()) { //< if not cached
        _taus[kpc] = apply<ParticleFinder>(currentEvent(), "Taus"+kpc).particlesByPt();
      }
      return _taus.at(kpc);
    }
    /// Get the pT-ordered set of reconstructed prompt (hadronic) taus, with extra cuts
    Particles taus(const Cut& cut, IDClass tauclass=IDClass::MEDIUM, bool include_nonprompt=false) const {
      return select(taus(tauclass, include_nonprompt), cut);
    }


    /// Get the pT-ordered set of tracks (charged particles)
    ///
    /// @note This less-used object type is not cached, as they are too low-level
    /// for e.g. overlap removal post-projection filtering.
    const Particles tracks() const {
      return apply<ParticleFinder>(currentEvent(), "Tracks").particlesByPt();
    }
    /// Get the pT-ordered set of tracks (charged particles), with extra cuts
    ///
    /// @note This less-used object type is not cached, as they are too low-level
    /// for e.g. overlap removal post-projection filtering.
    Particles tracks(const Cut& cut) const {
      return apply<ParticleFinder>(currentEvent(), "Tracks").particlesByPt(cut);
    }


    /// Get the set of jet-collection names (also usable for b-jets)
    const std::vector<std::string> jetnames() const {
      std::vector<std::string> rtn;
      rtn.reserve(_jets.size());
      for (const auto& kv : _jets) rtn.push_back(kv.first);
      return rtn;
    }

    /// Get a pT-ordered named set of reconstructed jets
    const Jets& jets(const std::string& jetname="DEFAULT") const {
      if (_jets.find(jetname) == _jets.end()) { //< not cached
        /// @todo Error handling if no projection or cut registered for that jet name
        _jets[jetname] = apply<JetFinder>(currentEvent(), "Jets"+jetname).jetsByPt(_jcuts.at(jetname));
      }
      return _jets.at(jetname);
    }
    /// Get a pT-ordered named set of reconstructed jets, with extra cuts
    Jets jets(const Cut& cut, const std::string& jetname="DEFAULT") const {
      return select(jets(jetname), cut);
    }


    /// Get a pT-ordered named set of reconstructed b-tagged jets
    const Jets& bjets(const std::string& jetname="DEFAULT") const {
      if (_bjets.find(jetname) == _bjets.end()) { //< not cached
        _bjets[jetname] = select(jets(jetname), hasBTag(_bcut, _bdeltaR));
      }
      return _bjets.at(jetname);
    }
    /// Get a pT-ordered named set of reconstructed b-tagged jets, with extra cuts
    Jets bjets(const Cut& cut, const std::string& jetname="DEFAULT") const {
      return select(bjets(jetname), cut);
    }


    /// Get the MET vector
    Vector3 vmet() const {
      return apply<METFinder>(currentEvent(), "MET").vectorMissingPt();
    }

    /// Get the MET scalar
    double met() const {
      return vmet().mod();
    }

    /// Get the SET scalar
    double set() const {
      return apply<METFinder>(currentEvent(), "MET").scalarEt();
    }

    /// Get the MET significance
    double metSignf() const {
      return apply<SmearedMET>(currentEvent(), "MET").missingEtSignf();
    }

    /// @}


    /// @name Event-manipulation functions
    /// @{

    /// @brief Perform a simple, in-place overlap removal between all physics objects
    ///
    /// @note This is not a specific experiment, physics group, or analysis' overlap
    /// removal scheme: in general, every analysis does something different and if you
    /// want an exact replication of the experimental procedure, you need to implement
    /// exactly (as far as possible at truth/smeared-level) the procedure used in the
    /// original analysis. This isn't Rivet-specific: it applies to any analysis preservation!
    ///
    /// This implementation will do a reasonable basic job of removing double-counting,
    /// though. It first removes any jets (including b-jets) within 0.2 of an electron or
    /// muon. Then removes photons too close to electrons. And finally removes any electron
    /// or muon within the radius of the remaining jet. Suggestions for improvement or
    /// flexible/user-friendly generalisation are welcome!
    ///
    /// @warning Note that this currently does no overlap-removal with respect to taus.
    ///
    /// @todo Anything specific to taus? (E.g. OR wrt jets based on #tracks?) To b-jets?
    void doSimpleOverlapRemoval(IDClass eclass=IDClass::MEDIUM, IDClass muclass=IDClass::MEDIUM,
                                IDClass phoclass=IDClass::MEDIUM, //IDClass tauclass=IDClass::MEDIUM,
                                bool include_nonprompt=false, const std::string& jetname="DEFAULT") {
      // Get non-const handles for filtering (and ensure objects exist)
      Jets& myjets = const_cast<Jets&>(jets(jetname));
      Jets& mybjets = const_cast<Jets&>(bjets(jetname));
      Particles& myes = const_cast<Particles&>(electrons(eclass, include_nonprompt));
      Particles& myms = const_cast<Particles&>(muons(muclass, include_nonprompt));
      //Particles& myts = const_cast<Particles&>(taus(tauclass, include_nonprompt));
      Particles& myys = const_cast<Particles&>(photons(phoclass, include_nonprompt));
      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(myjets, myes, 0.2);
      idiscardIfAnyDeltaRLess(myjets, myms, 0.2);
      idiscardIfAnyDeltaRLess(mybjets, myes, 0.2);
      idiscardIfAnyDeltaRLess(mybjets, myms, 0.2);
      // Remove all isolated photons within dR < 0.2 of a dressed electron
      /// @todo Also remove some electrons?
      idiscardIfAnyDeltaRLess(myys, myes, 0.2);
      // Remove dressed leptons within the radius of a remaining jet
      idiscardIfAnyDeltaRLess(myes, myjets, _jdefs.at(jetname).R);
      idiscardIfAnyDeltaRLess(myms, myjets, _jdefs.at(jetname).R);
      /// @todo Add anything for hadronic taus? Tau-jet overlap is probably "correct"
    }


    /// Re-cluster small-radius jets as large-radius ones, writing to a new jet collection
    // void doJetRecluster(/* jet collection names (in and out), recluster-alg enum and dR */) {
    //   /// @todo ...
    // }

    /// @}


    /// @name Stats helpers
    /// @{

    /// Helper function to elide the xsec/sumW scalefactor
    template<typename T>
    void scaleToIntLumi(T& ao, double intlumi) {
      scale(ao, intlumi*crossSection()/sumOfWeights());
    }

    /// @}


    /// @name Projection access
    /// @{

    /// Direct access to the underlying projection
    const ParticleFinder& electronsProj(IDClass eclass=IDClass::MEDIUM, bool include_nonprompt=false) {
      return getProjection<ParticleFinder>(_pckey(include_nonprompt, eclass, "Electrons"));
    }
    /// Direct access to the underlying projection
    const ParticleFinder& muonsProj(IDClass muclass=IDClass::MEDIUM, bool include_nonprompt=false) {
      return getProjection<ParticleFinder>(_pckey(include_nonprompt, muclass, "Muons"));
    }
    /// Direct access to the underlying projection
    const ParticleFinder& tausProj(IDClass tauclass=IDClass::MEDIUM, bool include_nonprompt=false) {
      return getProjection<ParticleFinder>(_pckey(include_nonprompt, tauclass, "Taus"));
    }
    /// Direct access to the underlying projection
    const ParticleFinder& photonsProj(IDClass phoclass=IDClass::MEDIUM, bool include_nonprompt=false) {
      return getProjection<ParticleFinder>(_pckey(include_nonprompt, phoclass, "Photons"));
    }
    /// Direct access to the underlying projection
    const ParticleFinder& tracksProj() {
      return getProjection<ParticleFinder>("Tracks");
    }
    /// Direct access to the underlying projection
    const JetFinder& jetsProj(const std::string& jetname="DEFAULT") {
      return getProjection<JetFinder>("Jets_"+jetname);
    }
    /// Direct access to the underlying projection
    const METFinder& metProj() { return getProjection<METFinder>("MET"); }

    /// @}


  private:

    /// @name Efficiency and smearing functions
    /// @{
    std::map<IDClass, ParticleEffFn> _effs_electron{ { IDClass::LOOSE, PARTICLE_EFF_ONE }, { IDClass::MEDIUM, PARTICLE_EFF_ONE }, { IDClass::TIGHT, PARTICLE_EFF_ONE } };
    std::map<IDClass, ParticleEffFn> _effs_muon{ { IDClass::LOOSE, PARTICLE_EFF_ONE }, { IDClass::MEDIUM, PARTICLE_EFF_ONE }, { IDClass::TIGHT, PARTICLE_EFF_ONE } };
    std::map<IDClass, ParticleEffFn> _effs_photon{ { IDClass::LOOSE, PARTICLE_EFF_ONE }, { IDClass::MEDIUM, PARTICLE_EFF_ONE }, { IDClass::TIGHT, PARTICLE_EFF_ONE } };
    std::map<IDClass, ParticleEffFn> _effs_tau{ { IDClass::LOOSE, PARTICLE_EFF_ONE }, { IDClass::MEDIUM, PARTICLE_EFF_ONE }, { IDClass::TIGHT, PARTICLE_EFF_ONE } };
    ParticleSmearFn _smear_electron{PARTICLE_SMEAR_IDENTITY}, _smear_muon{PARTICLE_SMEAR_IDENTITY}, _smear_photon{PARTICLE_SMEAR_IDENTITY}, _smear_tau{PARTICLE_SMEAR_IDENTITY};
    ParticleEffFn _eff_trk{PARTICLE_EFF_ONE};
    ParticleSmearFn _smear_trk{PARTICLE_SMEAR_IDENTITY};
    std::map<std::string, JetSmearFn> _smears_jet{ { "DEFAULT", JET_SMEAR_IDENTITY } };
    std::map<std::string, JetEffFn> _effs_btag{ { "DEFAULT", JET_BTAG_IDENTITY } };
    METSmearParamsFn _smearps_met{MET_SMEARPARAMS_IDENTITY};
    METSmearFn _smear_met{MET_SMEAR_IDENTITY};
    /// @}

    /// @name Kinematic requirements on physics objects
    /// @{
    Cut _acut = Cuts::abseta < 5.0; //< base acceptance, applied to all
    Cut _ecut = Cuts::pT > 10*GeV;
    Cut _mcut = Cuts::pT > 10*GeV;
    Cut _tcut = Cuts::pT > 10*GeV;
    Cut _trkcut = Cuts::abseta < 2.5;
    Cut _phocut = Cuts::pT > 10*GeV;
    std::map<std::string, Cut> _jcuts = {{ "DEFAULT",  Cuts::pT > 20*GeV && Cuts::abseta < 5.0 }};
    std::map<std::string, JetScheme> _jdefs = {{ "DEFAULT", JetScheme(JetAlg::ANTIKT, 0.4) }};
    /// @}


    /// @name Kinematic requirements on tagging particles
    ///
    /// @todo Make specific to jet name?
    /// @{
    Cut _bcut = Cuts::pT > 5*GeV && Cuts::abseta < 2.5;
    double _bdeltaR = 0.3;
    /// @}


    /// @name Physics-object collections
    ///
    /// @note Mutable to allow const-correct caching
    /// @{

    // Internal typedefs
    // using ParticlesMap = std::map<std::string, Particles>;
    // using JetsMap = std::map<std::string, Jets>;

    /// Named maps of prompt/inclusive particle collections
    mutable std::map<std::string, Particles> _electrons, _muons, _taus, _photons;

    /// Named maps of jet collections
    mutable std::map<std::string, Jets> _jets, _bjets;

    /// @}


    /// List of ID classes (and their string reprs) to iterate over
    //
    /// @note Following advice at https://qr.ae/p2EvHB
    const std::map<IDClass, std::string>& _idclasses() const {
      // static const std::vector<IDClass> _idclasses{ {IDClass::LOOSE, IDClass::MEDIUM, IDClass::TIGHT} };
      static const auto* map = new std::map<IDClass, std::string>{
        {IDClass::LOOSE, "Loose"}, {IDClass::MEDIUM, "Medium"}, {IDClass::TIGHT, "Tight"} };
      return *map;
    }
    const std::string& _idclassname(const IDClass& idc) const {
      return _idclasses().at(idc);
    }

    /// Create a promptness-IDclass string for keying maps
    std::string _pckey(bool include_nonprompt, IDClass idc, const std::string& prefix="") const {
      return prefix + (include_nonprompt ? "All" : "Prompt") + _idclassname(idc);
    }

  };


}


/// @def RIVET_DEFAULT_SIMPLEANALYSIS_CTOR
///
/// Preprocessor define to prettify the awkward derived constructor, with a name-string argument
#define RIVET_DEFAULT_SIMPLEANALYSIS_CTOR(clsname) clsname() : SimpleAnalysis(# clsname) {}


/// @def setDetSmearing
/// A macro for setting a whole family of eff/smearing functions at once
///
/// @todo Specify electron medium/tight, book both to allow control, or generalise
#define setDetSmearing(DET, BTAG) do {                                  \
  setElectronReco(ELECTRON_EFF_ ## DET ## _LOOSE, ELECTRON_SMEAR_ ## DET, IDClass::LOOSE) ; \
  setElectronReco(ELECTRON_EFF_ ## DET ## _MEDIUM, ELECTRON_SMEAR_ ## DET, IDClass::MEDIUM) ; \
  setElectronReco(ELECTRON_EFF_ ## DET ## _TIGHT, ELECTRON_SMEAR_ ## DET, IDClass::TIGHT) ; \
  setMuonReco(MUON_EFF_ ## DET ## _LOOSE, MUON_SMEAR_ ## DET, IDClass::LOOSE) ; \
  setMuonReco(MUON_EFF_ ## DET ## _MEDIUM, MUON_SMEAR_ ## DET, IDClass::MEDIUM) ; \
  setMuonReco(MUON_EFF_ ## DET ## _TIGHT, MUON_SMEAR_ ## DET, IDClass::TIGHT) ; \
  setPhotonReco(PHOTON_EFF_ ## DET ## _LOOSE, ELECTRON_SMEAR_ ## DET, IDClass::LOOSE) ; \
  setPhotonReco(PHOTON_EFF_ ## DET ## _MEDIUM, ELECTRON_SMEAR_ ## DET, IDClass::MEDIUM) ; \
  setPhotonReco(PHOTON_EFF_ ## DET ## _TIGHT, ELECTRON_SMEAR_ ## DET, IDClass::TIGHT) ; \
  setTauReco(TAU_EFF_ ## DET ## _LOOSE, TAU_SMEAR_ ## DET, IDClass::LOOSE) ; \
  setTauReco(TAU_EFF_ ## DET ## _MEDIUM, TAU_SMEAR_ ## DET, IDClass::MEDIUM) ; \
  setTauReco(TAU_EFF_ ## DET ## _TIGHT, TAU_SMEAR_ ## DET, IDClass::TIGHT) ; \
  setTrkReco(TRK_EFF_ ## DET, TRK_SMEAR_ ## DET) ;                      \
  for (const std::string& jname : jetnames())                           \
    setJetReco(JET_SMEAR_ ## DET, JET_BTAG_ ## DET ## _ ## BTAG, jname) ; \
  setMETReco(MET_SMEARPARAMS_ ## DET, MET_SMEAR_ ## DET) ;                                \
  } while (0)
