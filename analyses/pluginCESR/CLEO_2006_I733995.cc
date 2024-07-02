// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief anti-deuteron spectrum in upslion decays
  class CLEO_2006_I733995 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2006_I733995);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_p,1,1,1);
      book(_w_ups,"TMP/w_ups");
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& deut) {
      // deuteron id code
      static const long id = -1000010020;
      for (const Particle& p: mother.children()) {
        if (p.pid() == id) {
          deut.push_back(p);
        }
        else if (!p.children().empty()) {
          findDecayProducts(p, deut);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==553) ) {
        _w_ups->fill();
        Particles deut;
        findDecayProducts(ups, deut);
        if (deut.empty()) continue;
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV) {
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
        }
        for (const Particle& p : deut) {
          const double mom = boost.transform(p.mom()).p3().mod();
          _h_p->fill(mom);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = .817+.022;
      scale(_h_p, 1e5/fact / *_w_ups);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_p;
    CounterPtr _w_ups;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2006_I733995);

}
