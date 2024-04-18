#include "Rivet/Tools/ParticleName.hh"
#include "Rivet/Tools/Utils.hh"

namespace Rivet {
  namespace PID {


    // Initialise ParticleNames singleton pointer
    unique_ptr<ParticleNames> ParticleNames::_instance = nullptr;


    std::string ParticleNames::_particleName(PdgId pid) {
      if (_ids_names.find(pid) == _ids_names.end()) {
        //throw PidError("Particle ID '" + lexical_cast<string>(pid) + "' not known.");
        return lexical_cast<string>(pid);
      }
      return _ids_names[pid];
    }


    PdgId ParticleNames::_particleID(const std::string& pname) {
      if (_names_ids.find(pname) == _names_ids.end()) {
        std::string up = toUpper(pname);
        if (up == "P+" || up == "P") return PROTON;
        if (up == "P-" || up == "PBAR") return ANTIPROTON;
        if (up == "E-") return ELECTRON;
        if (up == "E+") return POSITRON;
        if (up == "GAMMA") return PHOTON;
        if (up == "N") return NEUTRON;
        if (up == "D") return DEUTERON;
        if (up == "AL") return ALUMINIUM;
        if (up == "CU") return COPPER;
        if (up == "XE") return XENON;
        if (up == "AU") return GOLD;
        if (up == "PB") return LEAD;
        if (up == "U") return URANIUM;
        try {
          PdgId rtn = lexical_cast<PdgId>(pname);
          return rtn;
        } catch (const bad_lexical_cast& blc) {
          //throw PidError("Particle name '" + pname + "' not known and could not be directly cast to a PDG ID.");
          return 0;
        }
      }
      return _names_ids[pname];
    }


  }
}
