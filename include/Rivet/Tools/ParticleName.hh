#ifndef RIVET_PARTICLENAME_HH
#define RIVET_PARTICLENAME_HH

#include "Rivet/Particle.fhh"
#include "Rivet/Tools/Exceptions.hh"

// cython includes termio.h on some systems,
// which #defines a B0 = 0 macro
#undef B0

namespace Rivet {

  namespace PID {


    /// Static const convenience particle ID names

    /// Special wildcard particle name
    constexpr PdgId ANY = 10000;

    /// @name Charged leptons
    /// @{
    constexpr PdgId ELECTRON = 11;
    constexpr PdgId POSITRON = -ELECTRON;
    constexpr PdgId EMINUS = ELECTRON;
    constexpr PdgId EPLUS = POSITRON;
    constexpr PdgId MUON = 13;
    constexpr PdgId ANTIMUON = -MUON;
    constexpr PdgId TAU = 15;
    constexpr PdgId ANTITAU = -TAU;
    /// @}

    /// @name Neutrinos
    /// @{
    constexpr PdgId NU_E = 12;
    constexpr PdgId NU_EBAR = -NU_E;
    constexpr PdgId NU_MU = 14;
    constexpr PdgId NU_MUBAR = -NU_MU;
    constexpr PdgId NU_TAU = 16;
    constexpr PdgId NU_TAUBAR = -NU_TAU;
    /// @}

    /// @name Bosons
    /// @{
    constexpr PdgId PHOTON = 22;
    constexpr PdgId GAMMA = PHOTON;
    constexpr PdgId GLUON = 21;
    constexpr PdgId WPLUSBOSON = 24;
    constexpr PdgId WMINUSBOSON = -WPLUSBOSON;
    constexpr PdgId WPLUS = WPLUSBOSON;
    constexpr PdgId WMINUS = WMINUSBOSON;
    constexpr PdgId WBOSON = WPLUSBOSON;
    constexpr PdgId Z0BOSON = 23;
    constexpr PdgId ZBOSON = Z0BOSON;
    constexpr PdgId Z0 = Z0BOSON;
    constexpr PdgId HIGGSBOSON = 25;
    constexpr PdgId HIGGS = HIGGSBOSON;
    constexpr PdgId H0BOSON = HIGGSBOSON;
    constexpr PdgId HBOSON = HIGGSBOSON;
    /// @}

    /// @name Quarks
    /// @{
    constexpr PdgId DQUARK = 1;
    constexpr PdgId UQUARK = 2;
    constexpr PdgId SQUARK = 3;
    constexpr PdgId CQUARK = 4;
    constexpr PdgId BQUARK = 5;
    constexpr PdgId TQUARK = 6;
    /// @}

    /// @name Nucleons
    /// @{
    constexpr PdgId PROTON = 2212;
    constexpr PdgId ANTIPROTON = -PROTON;
    constexpr PdgId PBAR = ANTIPROTON;
    constexpr PdgId NEUTRON = 2112;
    constexpr PdgId ANTINEUTRON = -NEUTRON;
    /// @}

    /// @name Light mesons
    /// @{
    constexpr PdgId PI0 = 111;
    constexpr PdgId PIPLUS = 211;
    constexpr PdgId PIMINUS = -PIPLUS;
    constexpr PdgId RHO0 = 113;
    constexpr PdgId RHOPLUS = 213;
    constexpr PdgId RHOMINUS = -RHOPLUS;
    constexpr PdgId K0L = 130;
    constexpr PdgId K0S = 310;
    constexpr PdgId K0 = 311;
    constexpr PdgId KPLUS = 321;
    constexpr PdgId KMINUS = -KPLUS;
    constexpr PdgId ETA = 221;
    constexpr PdgId ETAPRIME = 331;
    constexpr PdgId PHI = 333;
    constexpr PdgId OMEGA = 223;
    /// @}

    /// @name Charmonia
    /// @{
    constexpr PdgId ETAC = 441;
    constexpr PdgId JPSI = 443;
    constexpr PdgId PSI2S = 100443;
    /// @}

    /// @name Charm mesons
    /// @{
    constexpr PdgId D0 = 421;
    constexpr PdgId D0BAR = -421;
    constexpr PdgId DPLUS = 411;
    constexpr PdgId DMINUS = -DPLUS;
    constexpr PdgId DSTARPLUS = 413;
    constexpr PdgId DSTARMINUS = -DSTARPLUS;
    constexpr PdgId DSPLUS = 431;
    constexpr PdgId DSMINUS = -DSPLUS;
    /// @}

    /// @name Bottomonia
    /// @{
    constexpr PdgId ETAB = 551;
    constexpr PdgId UPSILON1S = 553;
    constexpr PdgId UPSILON2S = 100553;
    constexpr PdgId UPSILON3S = 200553;
    constexpr PdgId UPSILON4S = 300553;
    /// @}

    /// @name b mesons
    /// @{
    constexpr PdgId B0 = 511;
    constexpr PdgId B0BAR = -511;
    constexpr PdgId BPLUS = 521;
    constexpr PdgId BMINUS = -BPLUS;
    constexpr PdgId B0S = 531;
    constexpr PdgId BCPLUS = 541;
    constexpr PdgId BCMINUS = -BCPLUS;
    /// @}

    /// @name Baryons
    /// @{
    constexpr PdgId LAMBDA = 3122;
    constexpr PdgId SIGMA0 = 3212;
    constexpr PdgId SIGMAPLUS = 3222;
    constexpr PdgId SIGMAMINUS = 3112;
    constexpr PdgId SIGMAB = 5212;
    constexpr PdgId SIGMABPLUS = 5222;
    constexpr PdgId SIGMABMINUS = 5112;
    constexpr PdgId LAMBDACPLUS = 4122;
    constexpr PdgId LAMBDACMINUS = -4122;
    constexpr PdgId LAMBDAB = 5122;
    constexpr PdgId XI0 = 3322;
    constexpr PdgId XIMINUS = 3312;
    constexpr PdgId XIPLUS = -XIMINUS;
    constexpr PdgId XI0B = 5232;
    constexpr PdgId XIBMINUS = 5132;
    constexpr PdgId XI0C = 4132;
    constexpr PdgId XICPLUS = 4232;
    constexpr PdgId OMEGAMINUS = 3334;
    constexpr PdgId OMEGAPLUS = -OMEGAMINUS;
    constexpr PdgId OMEGABMINUS = 5332;
    constexpr PdgId OMEGA0C = 4332;
    /// @}

    /// @name Exotic/weird stuff
    /// @{
    constexpr PdgId REGGEON = 110;
    constexpr PdgId POMERON = 990;
    constexpr PdgId ODDERON = 9990;
    constexpr PdgId GRAVITON = 39;
    constexpr PdgId NEUTRALINO1 = 1000022;
    constexpr PdgId GRAVITINO = 1000039;
    constexpr PdgId GLUINO = 1000021;
    constexpr int BPRIME = 7;
    constexpr int TPRIME = 8;
    constexpr int LPRIME = 17;
    constexpr int NUPRIME = 18;
    // constexpr int DARKMATTERSCALAR = 1000051;
    // constexpr int DARKMATTERFERMION = 1000052;
    // constexpr int DARKMATTERVECTOR = 1000053;
    /// @todo Add axion, black hole remnant, etc. on demand
    /// @}

    /// @name Nuclei
    /// @{
    constexpr PdgId DEUTERON = 1000010020;
    constexpr PdgId ALUMINIUM = 1000130270;
    constexpr PdgId COPPER = 1000290630;
    constexpr PdgId XENON = 1000541290;
    constexpr PdgId GOLD = 1000791970;
    constexpr PdgId LEAD = 1000822080;
    constexpr PdgId URANIUM = 1000922380;
    /// @todo Add other nuclei on demand
    /// @}



    /// Handler for particle name code <-> string conversion
    ///
    /// @todo Hide this utility class or equiv in an unnamed namespace in a .cc file.
    class ParticleNames {
    public:

      static std::string particleName(PdgId pid) {
        if (!_instance) _instance = unique_ptr<ParticleNames>(new ParticleNames);
        return _instance->_particleName(pid);
      }

      static PdgId particleID(const std::string& pname) {
        if (!_instance) _instance = unique_ptr<ParticleNames>(new ParticleNames);
        return _instance->_particleID(pname);
      }

      std::string _particleName(PdgId pid);

      PdgId _particleID(const std::string& pname);


    private:

      ParticleNames() {
        // charged leptons
        _add_pid_name(ELECTRON, "ELECTRON");
        _add_pid_name(POSITRON, "POSITRON");
        _add_pid_name(MUON, "MUON");
        _add_pid_name(ANTIMUON, "ANTIMUON");
        _add_pid_name(TAU, "TAU");
        _add_pid_name(ANTITAU, "ANTITAU");
        // neutrinos
        _add_pid_name(NU_E, "NU_E");
        _add_pid_name(NU_EBAR, "NU_EBAR");
        _add_pid_name(NU_MU, "NU_MU");
        _add_pid_name(NU_MUBAR, "NU_MUBAR");
        _add_pid_name(NU_TAU, "NU_TAU");
        _add_pid_name(NU_TAUBAR, "NU_TAUBAR");
        // bosons
        _add_pid_name(GLUON, "GLUON");
        _add_pid_name(PHOTON, "PHOTON");
        _add_pid_name(WPLUSBOSON, "WPLUSBOSON");
        _add_pid_name(WMINUSBOSON, "WMINUSBOSON");
        _add_pid_name(ZBOSON, "ZBOSON");
        _add_pid_name(HIGGS, "HIGGS");
        // quarks
        _add_pid_name(DQUARK, "DOWN");
        _add_pid_name(UQUARK, "UP");
        _add_pid_name(SQUARK, "STRANGE");
        _add_pid_name(CQUARK, "CHARM");
        _add_pid_name(BQUARK, "BOTTOM");
        _add_pid_name(TQUARK, "TOP");
        _add_pid_name(-DQUARK, "ANTIDOWN");
        _add_pid_name(-UQUARK, "ANTIUP");
        _add_pid_name(-SQUARK, "ANTISTRANGE");
        _add_pid_name(-CQUARK, "ANTICHARM");
        _add_pid_name(-BQUARK, "ANTIBOTTOM");
        _add_pid_name(-TQUARK, "ANTITOP");
        // nucleons
        _add_pid_name(PROTON, "PROTON");
        _add_pid_name(ANTIPROTON, "ANTIPROTON");
        _add_pid_name(NEUTRON, "NEUTRON");
        _add_pid_name(ANTINEUTRON, "ANTINEUTRON");
        // light mesons
        _add_pid_name(PI0, "PI0");
        _add_pid_name(PIPLUS, "PIPLUS");
        _add_pid_name(PIMINUS, "PIMINUS");
        _add_pid_name(RHO0, "RHO0");
        _add_pid_name(RHOPLUS, "RHOPLUS");
        _add_pid_name(RHOMINUS, "RHOMINUS");
        _add_pid_name(K0, "K0");
        _add_pid_name(KPLUS, "KPLUS");
        _add_pid_name(KMINUS, "KMINUS");
        _add_pid_name(ETA, "ETA");
        _add_pid_name(PHI, "PHI");
        // charm mesons
        _add_pid_name(D0, "D0");
        _add_pid_name(D0BAR, "D0BAR");
        _add_pid_name(DPLUS, "DPLUS");
        _add_pid_name(DMINUS, "DMINUS");
        // B mesons
        _add_pid_name(B0, "B0");
        _add_pid_name(B0BAR, "B0BAR");
        _add_pid_name(BPLUS, "BPLUS");
        _add_pid_name(BMINUS, "BMINUS");
        // baryons
        _add_pid_name(LAMBDA, "LAMBDA");
        _add_pid_name(SIGMA0, "SIGMA0");
        _add_pid_name(SIGMAPLUS, "SIGMAPLUS");
        _add_pid_name(SIGMAMINUS, "SIGMAMINUS");
        _add_pid_name(XI0, "XI0");
        _add_pid_name(XIPLUS, "XIPLUS");
        _add_pid_name(XIMINUS, "XIMINUS");
        // exotics
        _add_pid_name(REGGEON, "REGGEON");
        _add_pid_name(POMERON, "POMERON");
        // nuclei
        _add_pid_name(DEUTERON, "DEUTERON");
        _add_pid_name(ALUMINIUM, "ALUMINIUM");
        _add_pid_name(COPPER, "COPPER");
        _add_pid_name(XENON, "XENON");
        _add_pid_name(GOLD, "GOLD");
        _add_pid_name(LEAD, "LEAD");
        _add_pid_name(URANIUM, "URANIUM");
        // wildcard
        _add_pid_name(ANY, "*");
      }

      void _add_pid_name(PdgId pid, const std::string& pname) {
        _ids_names[pid] = pname;
        _names_ids[pname] = pid;
      }

      static unique_ptr<ParticleNames> _instance;

      std::map<PdgId, std::string> _ids_names;

      std::map<std::string, PdgId> _names_ids;

    };


    /// Print a PdgId as a named string.
    inline std::string toParticleName(PdgId p) {
      return ParticleNames::particleName(p);
    }


    /// Print a PdgId as a named string.
    inline PdgId toParticleID(const std::string& pname) {
      return ParticleNames::particleID(pname);
    }


    /// Convenience maker of particle ID pairs from PdgIds.
    inline std::pair<PdgId,PdgId> make_pdgid_pair(PdgId a, PdgId b) {
      return make_pair(a, b);
    }


    /// Convenience maker of particle ID pairs from particle names.
    inline std::pair<PdgId,PdgId> make_pdgid_pair(const std::string& a, const std::string& b) {
      const PdgId pa = toParticleID(a);
      const PdgId pb = toParticleID(b);
      return make_pair(pa, pb);
    }


    /// Print a PdgIdPair as a string.
    inline std::string toBeamsString(const PdgIdPair& pair) {
      string out = "[" +
        toParticleName(pair.first) + ", " +
        toParticleName(pair.second) + "]";
      return out;
    }


  }

}

#endif
