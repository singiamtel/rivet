// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Upsilon(4s) -> B Bbar
  class ARGUS_1990_I296188 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1990_I296188);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histos
      book(_h, 1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
      	// boost to rest frame
      	LorentzTransform cms_boost;
      	// loop over decay products
      	for (const Particle& p2 : p.children()) {
      	  if (p2.abspid()==511 || p2.abspid()==521) {
	          _h->fill(abs(p2.p3().z()/p2.p3().mod()));
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1990_I296188);

}
