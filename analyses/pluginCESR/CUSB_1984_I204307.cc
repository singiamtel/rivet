// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Photon spectrum in Upsilon(1,2S) decays
  class CUSB_1984_I204307 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CUSB_1984_I204307);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==553 || Cuts::pid==100553), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
	      book(_h[ix],1+ix,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      for (const Particle& ups : apply<UnstableParticles>(event, "UFS").particles()) {
        unsigned int iups = ups.pid()/100000;
      	unsigned int nhadron(0);
      	Particles photons;
      	for (const Particle& child : ups.children()) {
      	  if (PID::isHadron(child.pid())) {
      	    ++nhadron;
          }
      	  else if(child.pid()==PID::PHOTON) {
      	    photons.push_back(child);
          }
      	}
        if (nhadron==0 || photons.empty()) continue;
        LorentzTransform boost;
        if (ups.p3().mod() > 1*MeV) {
          boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
        }
        for (const Particle& gamma : photons) {
          const double E = boost.transform(gamma.mom()).E();
          // E = randnorm(E,0.21*sqrt(E));
          _h[iups]->fill(2.*E/ups.mass());
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


  RIVET_DECLARE_PLUGIN(CUSB_1984_I204307);

}
