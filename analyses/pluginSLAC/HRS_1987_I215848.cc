// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Hadron Spectra in $e^+e^-$ collisions at 29 GeV
  class HRS_1987_I215848 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1987_I215848);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h["pi"],  2, 1, 1);
      book(_h["Kp"],  3, 1, 1);
      book(_h["p"],   4, 1, 1);
      book(_h["K0"],  5, 1, 1);
      book(_h["lam"], 6, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles()) {
        const double xE = 2.*p.E()/sqrtS();
        const double beta = p.p3().mod() / p.E();
        if (p.pid()==130 || p.pid()==310) {
          _h["K0"]->fill(xE,1./beta);
        }
        else if (p.abspid()==321) {
          _h["Kp"]->fill(xE, 1./beta);
        }
        else if (p.abspid()==211) {
          _h["pi"]->fill(xE, 1./beta);
        }
        else if (p.abspid()==2212) {
          _h["p"]->fill(xE, 1./beta);
        }
        else if (p.abspid()==3122) {
          _h["lam"]->fill(xE, 1./beta);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()*sqr(sqrtS())/microbarn/sumOfWeights();
      scale(_h, sf);
    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1987_I215848);

}
