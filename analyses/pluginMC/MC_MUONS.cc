// -*- C++ -*-
#include "Rivet/Analyses/MC_PARTICLES_BASE.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// @brief MC validation analysis for muons
  class MC_MUONS : public MC_PARTICLES_BASE {
  public:

    MC_MUONS()
      : MC_PARTICLES_BASE("MC_MUONS", 2, "muon")
    {    }


    void init() {
      const bool direct = getOption<bool>("DIRECT", false);
      const bool dressed = getOption<bool>("DRESSED", direct);
      MSG_DEBUG("Direct-only: " << direct << ", dressed: " << dressed);
      FinalState muons(Cuts::abspid == PID::MUON);
      if (!direct) {
        declare(muons, "Muons");
      } else if (!dressed) {
        declare(PromptFinalState(muons), "Muons");
      } else {
        LeptonFinder dleps(muons, FinalState(Cuts::abspid == PID::PHOTON), 0.1);
        declare(dleps, "Muons");
      }

      MC_PARTICLES_BASE::init();
    }


    void analyze(const Event& event) {
      const Particles mus = apply<ParticleFinder>(event, "Muons").particlesByPt(Cuts::pT > 0.5*GeV);
      MC_PARTICLES_BASE::_analyze(event, mus);
    }


    void finalize() {
      MC_PARTICLES_BASE::finalize();
    }

  };


  RIVET_DECLARE_PLUGIN(MC_MUONS);

}
