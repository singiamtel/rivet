// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class LHCB_2011_I919315 : public Analysis {
  public:
    /// @name Constructors etc.
    /// @{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2011_I919315);

    /// @}
  public:

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
    	// select phi mesons in fiducial phase-space (symmetrical LHCb)
      declare(UnstableParticles(Cuts::abspid == 333 && Cuts::absrapIn(2.44, 4.06) && Cuts::ptIn(0.6*GeV, 5.0*GeV)), "phiFS");

      book(_h_Phi_pT_y, {2.44, 2.62, 2.8, 2.98, 3.16, 3.34, 3.52, 3.7, 3.88, 4.06},
           { "d02-x01-y01", "d02-x01-y02", "d03-x01-y01", "d03-x01-y02",
             "d04-x01-y01", "d04-x01-y02", "d05-x01-y01", "d05-x01-y02", "d06-x01-y01"});
      book(_h_Phi_pT ,7, 1, 1);
      book(_h_Phi_y ,8, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze (const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles> (event, "phiFS");

      for (const Particle& p : ufs.particles()) {
      	double y  = p.absrapidity();
        double pT = p.pt();
        _h_Phi_y->fill(y);
        _h_Phi_pT->fill(pT/MeV);
        _h_Phi_pT_y->fill(y, pT/GeV);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
    	// correct for symmetrical detector
      double scale_factor = crossSectionPerEvent()/microbarn/2.;
      scale (_h_Phi_y, scale_factor);
      scale (_h_Phi_pT, scale_factor);
      scale(_h_Phi_pT_y, scale_factor/1000.);
      divByGroupWidth(_h_Phi_pT_y);
    }

    /// @}

  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_Phi_y;
    Histo1DPtr _h_Phi_pT;
    Histo1DGroupPtr _h_Phi_pT_y;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(LHCB_2011_I919315);

}

// @}
