// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D -> phi X
  class BESIII_2019_I1749793 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1749793);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==30443), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }

    void findChildren(const Particle & p, Particles & phi) {
      for (const Particle& child : p.children()) {
        if (child.pid()==333) {
          phi.push_back(child);
        }
        else if (!child.children().empty()) {
          findChildren(child,phi);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles psi = apply<UnstableParticles>(event, "UFS").particles();
      if (psi.size()!=1) vetoEvent;
      for (const Particle& D : psi[0].children()) {
        int ih=-1;
        if      (D.abspid()==411) ih=0;
        else if (D.abspid()==421) ih=1;
        else  continue;
        Particles phi;
        findChildren(D,phi);
        for (const Particle& p : phi) {
          _h[ih]->fill(p.momentum().p3().mod());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1749793);

}
