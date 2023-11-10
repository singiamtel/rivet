// -*- C++ -*-
#include "Rivet/Particle.hh"
#include "Rivet/DressedLepton.hh"

namespace Rivet {


  DressedLepton::DressedLepton(const Particle& dlepton)
    : Particle(dlepton)
  {
    if (dlepton.isComposite()) {
      Particles dressing;
      dressing.reserve(dlepton.constituents().size()-1);
      for (const Particle& p : dlepton.constituents()) {
        if (p.isChargedLepton()) setConstituents({{p}}); //< bare lepton is first constituent
        else dressing.push_back(p);
      }
      addConstituents(dressing);
    } else {
      setConstituents({{dlepton}});
    }
  }

  DressedLepton::DressedLepton(const Particle& lepton, const Particles& photons, bool momsum)
    : Particle(lepton.pid(), lepton.momentum())
  {
    setConstituents({{lepton}}); //< bare lepton is first constituent
    addConstituents(photons, momsum);
  }

  void DressedLepton::addPhoton(const Particle& p, bool momsum) {
    if (p.pid() != PID::PHOTON) throw Error("Clustering a non-photon on to a DressedLepton:"+to_string(p.pid()));
    addConstituent(p, momsum);
  }

  const Particle& DressedLepton::bareLepton() const {
    const Particle& l = constituents().front();
    if (!l.isChargedLepton()) throw Error("First constituent of a DressedLepton is not a bare lepton: oops");
    return l;
  }


}
