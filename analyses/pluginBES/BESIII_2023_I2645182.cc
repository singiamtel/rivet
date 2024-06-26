// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s -> f0 e+ nu_e
  class BESIII_2023_I2645182 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2645182);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==431), "UFS");
      // histos
      book(_h,1,1,1);
      book(_nDs, "TMP/nDs");
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int mesonID=9010221;
      // Loop over Ds mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        _nDs->fill();
        Particle f0;
        if (isSemileptonicDecay(p, {mesonID, PID::POSITRON, PID::NU_E})) {
          f0  = select(p.children(), Cuts::pid==mesonID)[0];
          if (f0.children().size()==2 &&
              f0.children()[1].pid() == -f0.children()[0].pid() &&
              f0.children()[0].abspid()==211) {
            _h->fill((p.mom()-f0.mom()).mass2());
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // D_s lifetime pdg2022 (ns)
      const double tau = 504e-6;
      scale(_h,1./tau/ *_nDs);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _nDs;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2645182);

}
